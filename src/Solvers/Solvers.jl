module Solvers

using LinearAlgebra

import Tulip.StandardForm


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

    # Centrality parameter
    μ::Tv
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

# Generic functions that all solvers must implement
# It is assumed that all pre-processing and pre-solve is done,
# and problem must be in standard form.
# Therefore, solvers should not modify the problem data.

#=
    A solver should have access to:
        - problem data (read-only) in standard form
        - working memory:
            - current solution
            - residuals
            - search direction
            - Cholesky factor / pre-conditionner
        - I/O channels (or a logger than handles it?)
=#


# Some basic interface functions, to be specialized for each solver type
"""
    get_termination_status(solver)

Return the solver's current status.
"""
function get_status(s::AbstractIPMSolver) end

"""
    get_solution_status(solver)

Return the primal solution status.
"""
function get_primal_solution_status(s::AbstractIPMSolver) end

"""
    get_dual_solution_status(solver)

Return the dual solution status.
"""
function get_dual_solution_status(s::AbstractIPMSolver) end


"""
    compute_residuals!(s, res, pt, A, b, c, uind, uval)

Compute residuals at given iterate.
"""
function compute_residuals!(
    s::AbstractIPMSolver,
    res::Residuals{Tv},
    pt::Point{Tv},
    A::AbstractMatrix{Tv}, b::Vector{Tv}, c::Vector{Tv},
    uind::Vector{Int}, uval::Vector{Tv}
) where{Tv<:Real}
    error("compute_residuals! not implemented for $(typeof(s)) solver.")
end

# Check stopping criteria
function check_stopping_criteria!(
    s::AbstractIPMSolver,

) end


include("HSDSolver/HSDSolver.jl")
# include("MPC/MPSSolver.jl")  # TODO

end  # module