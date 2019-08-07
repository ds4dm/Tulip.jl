"""
    HSDSolver

Solver for the homogeneous self-dual algorithm.
"""
mutable struct HSDSolver{Tv<:Real} <: AbstractIPMSolver{Tv}

    # =================
    #   Problem data
    # =================
    pb::StandardForm{Tv}


    # ================
    #   Book-keping
    # ================
    niter::Int  # Number of IPM iterations


    #=====================
        Working memory
    =====================#

    # ncon::Int    # Number of constraints
    # nvar::Int    # Number of variables
    # nvarub::Int  # Number of upper-bounded variables

    point::Point{Tv}    # Current primal-dual iterate
    res::Residuals{Tv}  # Residuals at current iterate

    # rxs::Tv  # rhs for complimentary products
    # rwz::Tv  # rhs for complimentary products
    # rtk::T          # rhs for homogeneous complimentary products

    # TODO: Constructor
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


"""
    compute_residuals!(::HSDSolver, res, pt, A, b, c, uind, uval)

In-place computation of primal-dual residuals at point `pt`.
"""
function compute_residuals!(
    ::HSDSolver,
    res::Residuals{Tv}, pt::Point{Tv},
    A::AbstractMatrix{Tv}, b::Vector{Tv}, c::Vector{Tv},
    uind::Vector{Int}, uval::Vector{Tv}
) where{Tv<:Real}
    # Primal residual
    # ``rp = t*b - A*x``
    mul!(res.rp, A, pt.x)
    rmul!(res.rp, -oneunit(Tv))
    axpy!(pt.t, b, res.rp)

    # Upper-bound residual
    # ``ru = t*u - w - x``
    rmul!(res.ru, zero(Tv))
    axpy!(-oneunit(Tv), pt.w, res.ru)
    @views axpy!(-oneunit(Tv), pt.x[uind], res.ru)
    axpy!(pt.t, uval, res.ru)

    # Dual residual
    # ``rd = t*c - A'*y - s + z``
    mul!(res.rd, transpose(A), pt.y)
    rmul!(res.rd, -oneunit(Tv))
    axpy!(pt.t, c, res.rd)
    axpy!(-oneunit(Tv), pt.s, res.rd)
    @views axpy!(oneunit(Tv), pt.z, res.rd[uind])

    # Gap residual
    res.rg = dot(c, pt.x) - dot(b, pt.y) - dot(uval, pt.z) + pt.k

    # Residuals norm
    res.rp_nrm = norm(res.rp, Inf)
    res.ru_nrm = norm(res.ru, Inf)
    res.rd_nrm = norm(res.rd, Inf)
    res.rg_nrm = norm(res.rg, Inf)

    return nothing
end



function check_stopping_criterion(hsd)