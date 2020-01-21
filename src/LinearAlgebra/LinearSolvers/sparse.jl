using SparseArrays
using SuiteSparse

# ==============================================================================
#   SparseIndefLinearSolver
# ==============================================================================

"""
    SparseIndefLinearSolver{Tv}

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses an LDLt factorization of the quasi-definite augmented system.
"""
mutable struct SparseIndefLinearSolver{Tv<:BlasReal} <: IndefLinearSolver{Tv}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Tv}
    θ::Vector{Tv}
    regP::Vector{Tv}  # primal regularization
    regD::Vector{Tv}  # dual regularization

    # Factorization
    F

    # TODO: constructor with initial memory allocation
    function SparseIndefLinearSolver(A::SparseMatrixCSC{Tv, Int}) where{Tv<:BlasReal}
        m, n = size(A)
        θ = ones(Tv, n)

        S = [
            spdiagm(0 => -θ)  A';
            A spdiagm(0 => ones(m))
        ]

        # TODO: PSD-ness checks
        F = ldlt(Symmetric(S))
        return new{Tv}(m, n, A, θ, ones(Tv, n), ones(Tv, m), F)
    end

end

"""
    update_linear_solver!(ls, θ, regP, regD)

Update LDLt factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update_linear_solver!(
    ls::SparseIndefLinearSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:BlasReal}
    # Sanity checks
    length(θ)  == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(regP) == ls.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(ls.n)"
    ))
    length(regD) == ls.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(ls.m)"
    ))

    # Update diagonal scaling
    ls.θ .= θ
    # Update regularizers
    ls.regP .= regP
    ls.regD .= regD

    # Re-compute factorization
    # TODO: Keep S in memory, only change diagonal
    S = [
        spdiagm(0 => -ls.θ .- regP)  ls.A';
        ls.A spdiagm(0 => regD)
    ]

    # TODO: PSD-ness checks
    ldlt!(ls.F, Symmetric(S))

    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::SparseIndefLinearSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:BlasReal}
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = ls.F \ ξ

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # Iterative refinement
    # TODO:
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = ls.A*dx + ls.regD .* dy - ξp
    # resD = - dx .* (ls.θ + ls.regP) + ls.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    # ξ1 = [resD; resP]
    # d1 = ls.F \ ξ1

    # # Update search direction
    # @views dx .-= d1[1:n]
    # @views dy .-= d1[(n+1):(m+n)]

    return nothing
end

# ==============================================================================
#   SparsePosDefLinearSolver
# ==============================================================================

"""
    SparsePosDefLinearSolver{Tv}

Linear solver for the 2x2 augmented system
```math
    [-(Θ^{-1} + Rp)   A'] [dx] = [xi_d]
    [   A             Rd] [dy] = [xi_p]
```
with ``A`` sparse.

Uses a Cholesky factorization of the positive definite normal equations system
```
(A*(Θ^{-1} + Rp)^{-1}*A' + Rd)  dy = xi_p + A*(Θ^{-1} + Rp)^{-1}*xi_d
                                dx = (Θ^{-1} + Rp)^{-1} * (A' dy - xi_d)
```
"""
mutable struct SparsePosDefLinearSolver{Tv<:BlasReal} <: PosDefLinearSolver{Tv}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Tv, Int}
    θ::Vector{Tv}   # Diagonal scaling
    # Regularization
    regP::Vector{Tv}  # primal
    regD::Vector{Tv}  # dual

    # Factorization
    F

    # Constructor and initial memory allocation
    # TODO: symbolic only + allocation
    function SparsePosDefLinearSolver(A::SparseMatrixCSC{Tv, Int}) where{Tv<:BlasReal}
        m, n = size(A)
        θ = ones(Tv, n)

        S = A * A' + spdiagm(0 => ones(m))

        # TODO: PSD-ness checks
        F = cholesky(Symmetric(S))

        return new{Tv}(m, n, A, θ, zeros(Tv, n), ones(Tv, m), F)
    end

end

"""
    update_linear_solver!(ls, θ, regP, regD)

Compute normal equation system matrix, and update the factorization.
"""
function update_linear_solver!(
    ls::SparsePosDefLinearSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:BlasReal}
    
    # Sanity checks
    length(θ) == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(regP) == ls.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(ls.n)"
    ))
    length(regD) == ls.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(ls.m)"
    ))

    ls.θ .= θ
    ls.regP .= regP  # Primal regularization is disabled for normal equations
    ls.regD .= regD

    # Re-compute factorization
    # D = (Θ^{-1} + Rp)^{-1}
    D = Diagonal(one(Tv) ./ (ls.θ .+ ls.regP))
    Rd = spdiagm(0 => ls.regD)
    S = ls.A * D * ls.A' + Rd

    # Update factorization
    cholesky!(ls.F, Symmetric(S), check=false)
    issuccess(ls.F) || throw(PosDefException(0))

    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::SparsePosDefLinearSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:BlasReal}
    m, n = ls.m, ls.n

    d = one(Tv) ./ (ls.θ .+ ls.regP)
    D = Diagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ ls.A * (D * ξd)

    # Solve augmented system
    dy .= (ls.F \ ξ_)

    # Recover dx
    dx .= D * (ls.A' * dy - ξd)

    # TODO: Iterative refinement
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = ls.A * dx + ls.regD .* dy - ξp
    # resD = - D \ dx + ls.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    
    return nothing
end