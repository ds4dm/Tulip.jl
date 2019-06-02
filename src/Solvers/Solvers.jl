module Solvers

using LinearAlgebra

"""
    AbstractSolver

Abstraction layer for IPM solvers.

An IPM solver implements an interior-point algorithm.
Currently supported:
    * Homogeneous self-dual (HSD)
    * Mehrotra predictor-corrector (MPC)
"""
abstract type AbstractSolver end

include("HSDSolver/HSDSolver.jl")
# include("MPC/MPSSolver.jl")  # TODO

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

Return the solver's termination status.
"""
function get_termination_status(s::AbstractSolver) end

"""
    get_solution_status(solver)

Return the primal solution status.
"""
function get_primal_solution_status(s::AbstractSolver) end

"""
    get_dual_solution_status(solver)

Return the dual solution status.
"""
function get_dual_solution_status(s::AbstractSolver) end




# Solve the optimization problem
# TODO: decide how much comes in and how much comes out
function optimize!(s::AbstractSolver) end

# Compute residuals for current iterate
# TODO: this function should be type-generic, to allow for dispatch
function compute_residuals!(s::AbstractSolver) end

# Check stopping criteria
function check_stopping_criteria!(s::AbstractSolver) end




end