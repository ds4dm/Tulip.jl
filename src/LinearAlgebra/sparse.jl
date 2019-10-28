using SparseArrays
using SuiteSparse

const SS_FLOAT = Union{Float32, Float64}

"""
    SparseLinearSolver{T}

"""
mutable struct SparseLinearSolver{T<:SS_FLOAT} <: TLPLinearSolver{T}
    # Up-to-date flag
    # Reset at the beginning of every IPM iteration
    up_to_date::Bool
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{T, Int}
    θ::Vector{T}

    # Factorization
    F

    # TODO: constructor with initial memory allocation
    function SparseLinearSolver(A::SparseMatrixCSC{T, Int}) where{T<:SS_FLOAT}
        
        m, n = size(A)
        θ = ones(T, n)

        S = [
            spdiagm(0 => -θ)  A';
            A spdiagm(0 => ones(m))
        ]

        # TODO: PSD-ness checks
        F = ldlt(Symmetric(S))

        return new{T}(false, m, n, A, θ, F)
    end

end


function update_linear_solver(ls::SparseLinearSolver{T}) where{T<:SS_FLOAT}
    ls.up_to_date && return

    # Re-compute factorization
    S = [
        spdiagm(0 => -(one(T) ./ ls.θ) .+ 1e-12)  ls.A';
        ls.A spdiagm(0 => 1e-12 .* ones(ls.m))
    ]

    # TODO: PSD-ness checks
    ldlt!(ls.F, Symmetric(S))

    ls.up_to_date = true
    
    return nothing
end


"""
    solve_augmented_system!(dx, dy, A, ls, θ, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
```
    -Θ*dx  + A'dy = ξd
     A*dx         = ξp
```

`A` and `θ` are given to perform iterative refinement if needed.

# Arguments
- `dx::Vector{Tv}, dy::Vector{Tv}`: Vectors of unknowns, modified in-place
- `ls::SparseLinearSolver{Tv}`: Linear solver for the augmented system
- `A::AbstractMatrix{Tv}`: Constraint matrix
- `θ::Vector{Tv}`: Diagonal scaling vector
- `ξp::Vector{Tv}, ξd::Vector{Tv}`: Right-hand-side vectors
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::SparseLinearSolver{Tv}, A::AbstractMatrix{Tv}, θ::Vector{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = size(A)
    
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
    rp = A*dx - ξp
    rd = -dx ./ θ + A' * dy - ξd

    ξ1 = [rd; rp]
    d1 = ls.F \ ξ1

    # Update search direction
    @views dx .-= d1[1:n]
    @views dy .-= d1[(n+1):(m+n)]

    return nothing
end