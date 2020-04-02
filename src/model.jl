mutable struct Solution{Tv}
    m::Int
    n::Int
    col_primal::Vector{Tv}
    col_dual::Vector{Tv}
    row_primal::Vector{Tv}
    row_dual::Vector{Tv}

    Solution{Tv}(m, n) where{Tv} = new{Tv}(
        m, n,
        Vector{Tv}(undef, n), Vector{Tv}(undef, n),
        Vector{Tv}(undef, m), Vector{Tv}(undef, m)
    )
end

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

    # IPM solver
    # If required, the problem is transformed to standard form
    # when instantiating the IPMSolver object.
    solver::Union{Nothing, AbstractIPMSolver{Tv}}

    # Problem solution
    solution::Union{Nothing, Solution{Tv}}

    # TODO: post-crushed solution (?)

    Model{Tv}() where{Tv} = new{Tv}(Parameters{Tv}(), ProblemData{Tv}(), nothing, nothing)
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

    # Instantiate solver
    model.solver = HSDSolver{Tv}(model.params, model.pbdata)
    # @info "HSD solver" model.solver.A model.solver.b model.solver.uind

    # Solve the problem
    if isnothing(model.solver)
        error("No IPM solver is attached to the model")
        return nothing
    else
        # TODO: add a try-catch for error handling
        optimize!(model.solver, model.params)
    end

    # Recover solution in original space
    sol = Solution{Tv}(model.pbdata.ncon, model.pbdata.nvar)
    _extract_solution!(sol, model.pbdata, model.solver)
    model.solution = sol

    # Done.
    return nothing
end

function _extract_solution!(sol::Solution{Tv}, pb::ProblemData{Tv}, hsd::HSDSolver{Tv}) where{Tv}

    # Extract column information
    τ = hsd.pt.t
    nfree = 0
    nvarupb = 0
    for (j, (l, u)) in enumerate(zip(pb.lvar, pb.uvar))
        # Recover primal variable and its reduced cost
        j_ = j + nfree
        if l == Tv(-Inf) && u == Tv(Inf)
            # free variable
            sol.col_primal[j] = (hsd.pt.x[j_] - hsd.pt.x[j_+1]) / τ
            sol.col_dual[j] = (hsd.pt.s[j_] - hsd.pt.s[j_+1]) / τ
            nfree += 1
        elseif l == Tv(-Inf) && isfinite(u)
            # Un-flip and push upper bound
            sol.col_primal[j] = u - hsd.pt.x[j_] / τ
            # Un-flip reduced cost
            sol.col_dual[j] = -hsd.pt.s[j_] / τ
        elseif isfinite(l) && isfinite(u)
            nvarupb += 1
            # Un-push lower bound
            sol.col_primal[j] = l + hsd.pt.x[j_] / τ
            # Reduced cost has two components
            sol.col_dual[j] = (hsd.pt.s[j_] - hsd.pt.z[nvarupb]) / τ
        else
            sol.col_primal[j] = l + hsd.pt.x[j_] / τ
            sol.col_dual[j] = hsd.pt.s[j_] / τ
        end
    end

    # Extract row information
    nslack = 0
    nslackupb = 0
    for i in 1:pb.ncon
        sol.row_dual[i] = hsd.pt.y[i] / τ
    end

    # Compute row primal
    for (i, row) in enumerate(pb.arows)
        ax = zero(Tv)
        for (j, aij) in zip(row.nzind, row.nzval)
            ax += aij * sol.col_primal[j]
        end
        sol.row_primal[i] = ax
    end

    return nothing
end