mutable struct Model{T}

    # Parameters
    params::Parameters{T}

    # TODO: model status
    #= Use an enum
        * Empty
        * Modified
        * OptimizationInProgress (optimize! is being called)
        * Solved (optimize! was called and the problem has not been modified since)
            TODO: some modifications should not change the solution status, e.g.:
                * changing names
                * changing objective constant
    =#
    status::TerminationStatus

    # Problem data
    pbdata::ProblemData{T}

    # Presolved problem
    # If presolved is disabled, this will point to m.pbdata
    presolve_data::Union{Nothing, PresolveData{T}}

    # IPM solver
    # If required, the problem is transformed to standard form
    # when instantiating the IPMSolver object.
    solver::Union{Nothing, AbstractIPMOptimizer{T}}

    # Problem solution (in original space)
    solution::Union{Nothing, Solution{T}}

    Model{T}() where{T} = new{T}(
        Parameters{T}(), Trm_NotCalled, ProblemData{T}(),
        nothing, nothing, nothing
    )
end

# TODO
# Basic functionalities (e.g., copy, empty, reset) should go in this file
# Interface-like should go in Interfaces
#=
    * optimize!
    * empty!
    * querying/setting parameters & attributes
    * build/modify problem through Model object
    * solution query
=#

import Base.empty!

function Base.empty!(m::Model{T}) where{T}
    m.pbdata = ProblemData{T}()
    m.status = Trm_NotCalled
    m.presolve_data = nothing
    m.solver = nothing
    m.solution = nothing

    return nothing
end

"""
    optimize!(model::Model{T})

Solve the optimization problem.
"""
function optimize!(model::Model{T}) where{T}

    # Set number of threads
    model.params.Threads >= 1 || error(
        "Number of threads must be > 0 (is $(model.params.Threads))"
    )
    BLAS.set_num_threads(model.params.Threads)

    # Print initial stats
    if model.params.OutputLevel > 0
        @printf "\nProblem info\n"
        @printf "  Name        : %s\n" model.pbdata.name
        @printf "  Constraints : %d\n" model.pbdata.ncon
        @printf "  Variables   : %d\n" model.pbdata.nvar
        @printf "  Non-zeros   : %d\n" sum(length.([col.nzind for col in model.pbdata.acols]))
    end

    pb_ = model.pbdata
    # Presolve
    # TODO: improve the if-else
    ps_options = model.params.Presolve
    if ps_options.Level > 0
        model.presolve_data = PresolveData(model.pbdata, ps_options)
        t_ = @elapsed st = presolve!(model.presolve_data)
        model.status = st

        if model.params.OutputLevel > 0
            ps = model.presolve_data
            nz0 = mapreduce(col -> length(col.nzind), +, model.pbdata.acols)
            nz = sum(ps.nzrow[ps.rowflag])
            @printf "\nReduced problem info\n"
            @printf "  Constraints : %d  (removed %d)\n" ps.nrow (ps.pb0.ncon - ps.nrow)
            @printf "  Variables   : %d  (removed %d)\n" ps.ncol (ps.pb0.nvar - ps.ncol)
            @printf "  Non-zeros   : %d  (removed %d)\n" nz (nz0 - nz)
            @printf "Presolve time : %.3fs\n" t_
        end

        # Check presolve status
        if st == Trm_Optimal || st == Trm_PrimalInfeasible || st == Trm_DualInfeasible || st == Trm_PrimalDualInfeasible
            model.params.OutputLevel > 0 && println("Presolve solved the problem.")

            # Perform post-solve
            sol0 = Solution{T}(model.pbdata.ncon, model.pbdata.nvar)
            postsolve!(sol0, model.presolve_data.solution, model.presolve_data)
            model.solution = sol0

            # Book-keeping
            # TODO: have a ModelStatus that indicates that model was solved by presolve
            model.solver = nothing

            # Done.
            return
        end

        # Presolve was not able to solve the problem
        extract_reduced_problem!(model.presolve_data)
        pb_ = model.presolve_data.pb_red
    end

    # Extract data in IPM form
    dat = IPMData(pb_, model.params.MatrixFactory)

    # Instantiate the IPM solver
    model.params.IPM.OutputLevel = model.params.OutputLevel
    model.solver = instantiate(model.params.IPM.Factory, dat, model.params.KKT)

    # Solve the problem
    # TODO: add a try-catch for error handling
    ipm_optimize!(model.solver, model.params.IPM)

    # Recover solution in original space
    sol_inner = Solution{T}(pb_.ncon, pb_.nvar)
    _extract_solution!(sol_inner, pb_, model.solver)

    # Post-solve
    if ps_options.Level > 0
        sol_outer = Solution{T}(model.pbdata.ncon, model.pbdata.nvar)
        postsolve!(sol_outer, sol_inner, model.presolve_data)
        model.solution = sol_outer
    else
        model.solution = sol_inner
    end

    model.status = model.solver.solver_status

    # Done.
    return nothing
end

function _extract_solution!(sol::Solution{T},
    pb::ProblemData{T}, ipm::AbstractIPMOptimizer{T}
) where{T}

    m, n = pb.ncon, pb.nvar

    # Extract column information
    # TODO: check for ray vs vertex
    sol.primal_status = ipm.primal_status
    sol.dual_status = ipm.dual_status

    is_primal_ray = (sol.primal_status == Sln_InfeasibilityCertificate)
    is_dual_ray = (sol.dual_status == Sln_InfeasibilityCertificate)
    sol.is_primal_ray = is_primal_ray
    sol.is_dual_ray = is_dual_ray
    τ_ = (is_primal_ray || is_dual_ray) ? one(T) : inv(ipm.pt.τ)

    @. sol.x = ipm.pt.x[1:n] * τ_
    @. sol.s_lower = ipm.pt.zl[1:n] * τ_
    @. sol.s_upper = ipm.pt.zu[1:n] * τ_

    # Extract row information
    @. sol.y_lower = pos_part.(ipm.pt.y) * τ_
    @. sol.y_upper = neg_part.(ipm.pt.y) * τ_

    # Compute row primal
    for (i, row) in enumerate(pb.arows)
        ax = zero(T)
        for (j, aij) in zip(row.nzind, row.nzval)
            ax += aij * sol.x[j]
        end
        sol.Ax[i] = ax
    end

    # Primal and dual objectives
    if sol.primal_status == Sln_InfeasibilityCertificate
        # Unbounded ray
        sol.z_primal = -T(Inf)
        sol.z_dual = -T(Inf)
    elseif sol.primal_status == Sln_Optimal || sol.primal_status == Sln_FeasiblePoint
        sol.z_primal = ipm.primal_objective
    else
        # Unknown solution status
        sol.z_primal = NaN
    end

    if sol.dual_status == Sln_InfeasibilityCertificate
        # Farkas proof of infeasibility
        sol.z_primal = T(Inf)
        sol.z_dual = T(Inf)
    elseif sol.dual_status == Sln_Optimal || sol.dual_status == Sln_FeasiblePoint
        # Dual solution is feasible
        sol.z_dual = ipm.dual_objective
    else
        # Unknown solution status
        sol.z_dual = NaN
    end

    return nothing
end
