using Printf

include("ipmdata.jl")

"""
    Point{T, Tv}

Primal-dual point.
"""
mutable struct Point{T, Tv}
    # Dimensions
    m::Int  # Number of constraints
    n::Int  # Number of variables
    p::Int  # Total number of finite variable bounds (lower and upper)

    # Primal variables
     x::Tv  # Original variables
    xl::Tv  # Lower-bound slack: `x - xl == l` (zero if `l == -∞`)
    xu::Tv  # Upper-bound slack: `x + xu == u` (zero if `u == +∞`)

    # Dual variables
     y::Tv  # Dual variables
    zl::Tv  # Lower-bound dual, zero if `l == -∞`
    zu::Tv  # Upper-bound dual, zero if `u == +∞`

    # HSD variables, only used with homogeneous form
    # Otherwise, one must ensure that (τ, κ) = (1, 0)
    hflag::Bool  # Is homogeneous embedding used?
    τ::T
    κ::T

    # Centrality parameter
    μ::T

    # Constructor
    Point{T, Tv}(m, n, p; hflag::Bool) where{T, Tv<:AbstractVector{T}} = new{T, Tv}(
        m, n, p,
        tzeros(Tv, n), tzeros(Tv, n), tzeros(Tv, n),
        tzeros(Tv, m), tzeros(Tv, n), tzeros(Tv, n),
        hflag, one(T), one(T),
        one(T)
    )
end

function update_mu!(pt::Point)
    pt.μ = (dot(pt.xl, pt.zl) + dot(pt.xu, pt.zu) + pt.hflag * (pt.τ * pt.κ)) / (pt.p + pt.hflag)
    return nothing
end


"""
    Residuals{T, Tv}

Data structure for IPM residual vectors.
"""
mutable struct Residuals{T, Tv}
    # Primal residuals
    rp::Tv  # rp = τ*b - A*x
    rl::Tv  # rl = τ*l - (x - xl)
    ru::Tv  # ru = τ*u - (x + xu)

    # Dual residuals
    rd::Tv  # rd = τ*c - (A'y + zl - zu)
    rg::T  # rg = c'x - (b'y + l'zl - u'zu) + κ

    # Residuals' norms
    rp_nrm::T  # |rp|
    rl_nrm::T  # |rl|
    ru_nrm::T  # |ru|
    rd_nrm::T  # |rd|
    rg_nrm::T  # |rg|
end


"""
    AbstractIPMOptimizer

Abstraction layer for IPM solvers.

An IPM solver implements an interior-point algorithm.
Currently supported:
    * Homogeneous self-dual (HSD)
"""
abstract type AbstractIPMOptimizer{T} end

"""
    ipm_optimize!(ipm)

Run the interior-point optimizer of `ipm`.
"""
function ipm_optimize! end

include("HSD/HSD.jl")
include("MPC/MPC.jl")
