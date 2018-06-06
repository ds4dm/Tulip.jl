import MathProgBase
const MPB = MathProgBase


export TulipSolver

"""
    TulipSolver

"""
struct TulipSolver <: MPB.AbstractMathProgSolver
    #=======================================================
        Optimization-related parameters
    =======================================================#

    # I/O behaviour
    output_level::Int       # 0 means no output, 1 means normal

    # Termination criteria
    n_iter_max::Int         # Maximum number of barrier iterations
    time_limit::Float64     # Time limit (in seconds)
    eps_tol_p::Float64      # Numerical tolerance for primal feasibility
    eps_tol_d::Float64      # Numerical tolerance for dual feasibility
    eps_tol_g::Float64      # Numerical tolerance for optimality gap
end

function TulipSolver(;
    output_level=1,
    n_iter_max=100,
    time_limit=Inf,
    eps_tol_p=10.0^-8,
    eps_tol_d=10.0^-8,
    eps_tol_g=10.0^-8
)

    TulipSolver(output_level, n_iter_max, time_limit, eps_tol_p, eps_tol_d, eps_tol_g)
end

"""
    TulipMathProgModel

Wrapper for MathProgBase models.
"""
mutable struct TulipMathProgModel <: MPB.AbstractLinearQuadraticModel
    inner::Model
end


function MPB.LinearQuadraticModel(s::TulipSolver)
    # create an empty model

    inner_model = Tulip.Model(
        output_level=s.output_level,
        n_iter_max=s.n_iter_max,
        time_limit=s.time_limit,
        eps_tol_p=s.eps_tol_p,
        eps_tol_d=s.eps_tol_d,
        eps_tol_g=s.eps_tol_g
    )
    
    TulipMathProgModel(inner_model)
end

function MPB.optimize!(m::TulipMathProgModel)

    optimize!(m.inner)

    return m.inner.status

end