mutable struct Model{Ta<:AbstractMatrix{<:Real}, Tv<:AbstractVector{<:Real}, Ti<:AbstractVector{<:Integer}}
    A::Ta
    b::Tv
    c::Tv
    uind::Ti
    uval::Tv
end

"""
    HSDSolver

Solver for the homogeneous self-dual algorithm.

Solves problems of the form
    ``
    \min_{x} c'x
    s.t.     A x = b
               x ⩽ u
               x ⩾ 0
    ``
"""
mutable struct HSDSolver{T<:Real, Tv<:AbstractVector{T}} <: AbstractSolver


    #=====================
        Working memory
    =====================#

    # ncon::Int    # Number of constraints
    # nvar::Int    # Number of variables
    # nvarub::Int  # Number of upper-bounded variables

    # Primal-dual iterate
    sol::Point{T, Tv}       # Current primal-dual solution
    best_sol::Point{T, Tv}  # Best (feasible) solution seen so far

    # Residuals
    rp::Tv      # Primal residual
    rp_norm::T  # Infinite norm of primal residual
    rd::Tv      # Dual residual
    rd_norm::T  # Infinite norm of dual residul
    ru::Tv      # Primal upper-bound residual
    ru_norm::T  # Infinite norm of primal upper-bound residual
    rg::T       # Gap residual

    rxs::Tv  # rhs for complimentary products
    rwz::Tv  # rhs for complimentary products
    rtk::T          # rhs for homogeneous complimentary products

end

function optimize!(
    solver::HSDSolver
)
    # If none provided, compute initial point


    # Compute residuals

    # Log for iteration 0

    # Check for stopping criteria

    # Main loop
    keep_going = true
    while keep_going

        keep_going = false

    end

    # Termination status

    # Solution status

    # Return

    return true

end

function compute_residuals!(model::Model, solver::HSDSolver)

    # Primal residual
    # ``rp = t*b - A*x``
    mul!(solver.rp, model.A, solver.x)
    rmul!(solver.rp, -oneunit(eltype(solver.rp)))
    axpy!(solver.t, model.b, solver.rp)

    # Upper-bound residual
    # ``ru = t*u - w - x``
    rmul!(solver.ru, -zero(eltype(solver.ru)))
    axpy!(-oneunit(eltype(solver.x)), solver.w, solver.ru)
    @views axpy!(-oneunit(eltype(solver.x)), solver.x[model.uind], solver.ru)
    axpy!(solver.t, model.uval, solver.ru)

    # Dual residual
    # ``rd = t*c - A'*y - s + z``
    mul!(solver.rd, transpose(model.A), solver.y)
    rmul!(solver.rd, -oneunit(eltype(solver.rd)))
    axpy!(solver.t, model.c, solver.rd)
    axpy!(-oneunit(eltype(solver.s)), solver.s, solver.rd)
    @views axpy!(oneunit(eltype(solver.z)), solver.z, solver.ru[model.uind])

    # Gap residual
    solver.rg = dot(model.c, solver.x) - dot(model.b, solver.y) - dot(model.uval, solver.z) + solver.k

    return nothing
end