using Printf

"""
    Point{Tv<:Real}

Primal-dual point.
"""
mutable struct Point{Tv<:Real}
    # Dimensions
    m::Int  # Number of constraints
    n::Int  # Number of variables
    p::Int  # Number of upper-bounded variables

    # Primal variables
    x::Vector{Tv}  # Original variables
    w::Vector{Tv}  # Upper-bound slacks
    t::Tv          # Homogeneous variable

    # Dual variables
    y::Vector{Tv}  # Dual variables
    s::Vector{Tv}  # Reduced costs
    z::Vector{Tv}  # Dual of upper-bound
    k::Tv          # Homogeneous variable

    # Regularizers
    qp::Vector{Tv}  # Primal regularizer
    qd::Vector{Tv}  # Dual regularizer

    # Centrality parameter
    μ::Tv

    Point{Tv}(ncon::Int, nvar::Int, nupb::Int) where{Tv<:Real} = new{Tv}(
        ncon, nvar, nupb,
        zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zeros(Tv, ncon), zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zeros(Tv, nvar), zeros(Tv, ncon),
        zero(Tv)
    )
end

function update_mu!(pt::Point{Tv}) where{Tv<:Real}
    pt.μ = ((dot(pt.x, pt.s) + dot(pt.w, pt.z) + pt.t * pt.k)) / (pt.n + pt.p + 1)
    return nothing
end


"""
    Residuals{Tv<:Real}

"""
mutable struct Residuals{Tv<:Real}
    # Primal residuals
    rp::Vector{Tv}  # rp = t*b - A*x
    ru::Vector{Tv}  # ru = t*u - w - U*x

    # Dual residuals
    rd::Vector{Tv}  # rd = t*c - A'y - s + U'z
    rg::Tv          # rg = c'x - b'y - u'z + k

    # Residuals' norms
    rp_nrm::Tv  # |rp|
    ru_nrm::Tv  # |ru|
    rd_nrm::Tv  # |rd|
    rg_nrm::Tv  # |rg|
end


"""
    AbstractIPMSolver

Abstraction layer for IPM solvers.

An IPM solver implements an interior-point algorithm.
Currently supported:
    * Homogeneous self-dual (HSD)
    * Mehrotra predictor-corrector (MPC)
"""
abstract type AbstractIPMSolver{Tv<:Real} end

include("HSDSolver/HSDSolver.jl")
# TODO: MPC solver