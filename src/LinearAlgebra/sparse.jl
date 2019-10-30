using SparseArrays
using SuiteSparse

"""
    SparseLinearSolver{T}

"""
mutable struct SparseLinearSolver{T<:BlasReal} <: TLPLinearSolver{T}
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
    function SparseLinearSolver(A::SparseMatrixCSC{T, Int}) where{T<:BlasReal}
        m, n = size(A)
        θ = ones(T, n)

        S = [
            spdiagm(0 => -θ)  A';
            A spdiagm(0 => ones(m))
        ]

        # TODO: PSD-ness checks
        F = ldlt(Symmetric(S))

        return new{T}(m, n, A, θ, F)
    end

end

TLPLinearSolver(A::SparseMatrixCSC{Tv, Int64}) where{Tv<:BlasReal} = SparseLinearSolver(A)

function update_linear_solver(
    ls::SparseLinearSolver{Tv},
    d::AbstractVector{Tv}
) where{Tv<:BlasReal}
    
    # Sanity checks
    length(d) == ls.n || throw(DimensionMismatch(
        "d has length $(length(d)) but linear solver is for n=$(ls.n)."
    ))

    ls.θ .= d

    # Re-compute factorization
    S = [
        spdiagm(0 => -(one(Tv) ./ ls.θ) .+ 1e-6)  ls.A';
        ls.A spdiagm(0 => 1e-6 .* ones(ls.m))
    ]

    # TODO: PSD-ness checks
    ldlt!(ls.F, Symmetric(S))
    
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
    ls::SparseLinearSolver{Tv}, A::SparseMatrixCSC{Tv}, θ::Vector{Tv},
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