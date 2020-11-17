"""
    Point{T, Tv}

Primal-dual point.
"""
mutable struct Point{T, Tv}
    # Dimensions
    m::Int  # Number of constraints
    n::Int  # Number of variables
    p::Int  # Total number of finite variable bounds (lower and upper)
    hflag::Bool  # Is homogeneous embedding used?

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
    τ::T
    κ::T

    # Centrality parameter
    μ::T

    # Constructor
    Point{T, Tv}(m, n, p; hflag::Bool) where{T, Tv<:AbstractVector{T}} = new{T, Tv}(
        m, n, p, hflag,
        # Primal variables
        tzeros(Tv, n), tzeros(Tv, n), tzeros(Tv, n),
        # Dual variables
        tzeros(Tv, m), tzeros(Tv, n), tzeros(Tv, n),
        # Homogeneous variables
        one(T), one(T),
        # Centrality parameter
        one(T)
    )
end

function update_mu!(pt::Point)
    pt.μ = (dot(pt.xl, pt.zl) + dot(pt.xu, pt.zu) + pt.hflag * (pt.τ * pt.κ)) / (pt.p + pt.hflag)
    return nothing
end