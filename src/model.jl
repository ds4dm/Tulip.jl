mutable struct Model{Tv}

    # Parameters
    params::Parameters{Tv}

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

    # Problem data
    pbdata::ProblemData{Tv}

    # Presolved problem
    # If presolved is disabled, this will point to m.pbdata
    presolve_data::Union{Nothing, PresolveData{Tv}}

    # IPM solver
    # If required, the problem is transformed to standard form
    # when instantiating the IPMSolver object.
    solver::Union{Nothing, AbstractIPMSolver{Tv}}

    # Problem solution (in original space)
    solution::Union{Nothing, Solution{Tv}}

    Model{Tv}() where{Tv} = new{Tv}(
        Parameters{Tv}(), ProblemData{Tv}(),
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

function Base.empty!(m::Model{Tv}) where{Tv}
    m.pbdata = ProblemData{Tv}()
    m.presolve_data = nothing
    m.solver = nothing
    m.solution = nothing
    
    return nothing
end

"""
    optimize!(model::Model{Tv})

Solve the optimization problem.
"""
function optimize!(model::Model{Tv}) where{Tv}
    # @info model.pbdata
    
    # Presolve
    # TODO: presolve flag
    model.presolve_data = PresolveData(model.pbdata)
    st = presolve!(model.presolve_data)

    # Check presolve status
    st = model.presolve_data.status
    if st == Trm_Optimal || st == Trm_PrimalInfeasible || st == Trm_DualInfeasible || st == Trm_PrimalDualInfeasible
        @info "Presolve solved the problem"
        
        # Perform post-solve
        sol0 = Solution{Tv}(model.pbdata.ncon, model.pbdata.nvar)
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
    presolved_problem = model.presolve_data.pb_red

    # Instantiate solver
    model.solver = HSDSolver{Tv}(model.params, presolved_problem)

    # Solve the problem
    # TODO: add a try-catch for error handling
    optimize!(model.solver, model.params)

    # TODO: check solver's status

    # Recover solution in presolve space
    sol_inner = Solution{Tv}(presolved_problem.ncon, presolved_problem.nvar)
    _extract_solution!(sol_inner, presolved_problem, model.solver)

    # Post-solve
    sol_outer = Solution{Tv}(model.pbdata.ncon, model.pbdata.nvar)
    postsolve!(sol_outer, sol_inner, model.presolve_data)
    model.solution = sol_outer

    # Done.
    return nothing
end

function _extract_solution!(sol::Solution{Tv}, pb::ProblemData{Tv}, hsd::HSDSolver{Tv}) where{Tv}

    # Extract column information
    # TODO: check for ray vs vertex
    sol.primal_status = hsd.primal_status
    sol.dual_status = hsd.dual_status

    is_primal_ray = (sol.primal_status == Sln_InfeasibilityCertificate)
    is_dual_ray = (sol.dual_status == Sln_InfeasibilityCertificate)
    τ_ = (is_primal_ray || is_dual_ray) ? one(Tv) : inv(hsd.pt.t)
    
    nfree = 0
    nvarupb = 0
    for (j, (l, u)) in enumerate(zip(pb.lvar, pb.uvar))
        # Recover primal variable and its reduced cost
        j_ = j + nfree
        if l == Tv(-Inf) && u == Tv(Inf)
            # free variable
            sol.x[j] = (hsd.pt.x[j_] - hsd.pt.x[j_+1]) * τ_
            s = (hsd.pt.s[j_] - hsd.pt.s[j_+1]) * τ_
            nfree += 1
        elseif l == Tv(-Inf) && isfinite(u)
            # Un-flip and push upper bound
            sol.x[j] = (!is_primal_ray * u) - hsd.pt.x[j_] * τ_
            # Un-flip reduced cost
            s = -hsd.pt.s[j_] * τ_
        elseif isfinite(l) && isfinite(u)
            nvarupb += 1
            # Un-push lower bound
            sol.x[j] = (!is_primal_ray * l) + hsd.pt.x[j_] * τ_
            # Reduced cost has two components
            s = (hsd.pt.s[j_] - hsd.pt.z[nvarupb]) * τ_
        else
            sol.x[j] = (!is_primal_ray * l) + hsd.pt.x[j_] * τ_
            s = hsd.pt.s[j_] * τ_
        end

        # Reduced cost
        sol.s_lower[j] = pos_part(s)
        sol.s_upper[j] = neg_part(s)
    end

    # Extract row information
    nslack = 0
    nslackupb = 0
    for i in 1:pb.ncon
        y = hsd.pt.y[i] * τ_
        sol.y_lower[i] = pos_part(y)
        sol.y_upper[i] = neg_part(y)
    end

    # Compute row primal
    for (i, row) in enumerate(pb.arows)
        ax = zero(Tv)
        for (j, aij) in zip(row.nzind, row.nzval)
            ax += aij * sol.x[j]
        end
        sol.Ax[i] = ax
    end

    # Primal and dual objectives
    if sol.primal_status == Sln_InfeasibilityCertificate
        # Unbounded ray
        sol.z_primal = -Inf
        sol.z_dual = -Inf
    elseif sol.primal_status == Sln_Optimal || sol.primal_status == Sln_FeasiblePoint
        sol.z_primal = hsd.primal_bound_scaled
    else
        # Unknown solution status
        sol.z_primal = NaN
    end

    if sol.dual_status == Sln_InfeasibilityCertificate
        # Farkas proof of infeasibility
        sol.z_primal = Inf
        sol.z_dual = Inf
    elseif sol.dual_status == Sln_Optimal || sol.dual_status == Sln_FeasiblePoint
        # Dual solution is feasible
        sol.z_dual = hsd.dual_bound_scaled

        # sol.z_dual = pb.obj0
        # for (yl, yu, l, u) in zip(sol.y_lower, sol.y_upper, pb.lcon, pb.ucon)
        #     isfinite(l) && (sol.z_dual += yl * l)
        #     isfinite(u) && (sol.z_dual -= yu * u)
        # end
        # for (sl, su, l, u) in zip(sol.s_lower, sol.s_upper, pb.lvar, pb.uvar)
        #     isfinite(l) && (sol.z_dual += sl * l)
        #     isfinite(u) && (sol.z_dual -= su * u)
        # end
    else
        # Unknown solution status
        sol.z_dual = NaN
    end

    return nothing
end