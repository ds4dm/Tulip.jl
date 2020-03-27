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

    # TODO: post-crushed solution (?)

    Model{Tv}() where{Tv} = new{Tv}(Parameters{Tv}(), ProblemData{Tv}(), nothing)
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
    
    return nothing
end