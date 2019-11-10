using SparseArrays
using SuiteSparse

construct_matrix(
    ::Type{SparseMatrixCSC}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{Tv}
) where{Tv<:Real} = sparse(aI, aJ, aV, m, n)

# ==============================================================================
#   SparseIndefLinearSolver
# ==============================================================================

"""
    SparseIndefLinearSolver{Tv}

"""
mutable struct SparseIndefLinearSolver{Tv<:BlasReal, Ta<:SparseMatrixCSC{Tv, <:Integer}} <: IndefLinearSolver{Tv, Ta}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::Ta
    θ::Vector{Tv}
    rp::Vector{Tv}  # primal regularization
    rd::Vector{Tv}  # dual regularization

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
        return new{Tv, SparseMatrixCSC{Tv, Int}}(m, n, A, θ, ones(Tv, n), ones(Tv, m), F)
    end

end

AbstractLinearSolver(A::SparseMatrixCSC{Tv, Int64}) where{Tv<:BlasReal} = SparsePosDefLinearSolver(A)

function update_linear_solver!(
    ls::SparseIndefLinearSolver{Tv},
    θ::AbstractVector{Tv},
    rp::AbstractVector{Tv}=zeros(Tv, ls.n),
    rd::AbstractVector{Tv}=zeros(Tv, ls.m)
) where{Tv<:BlasReal}
    # Sanity checks
    length(θ)  == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(rp) == ls.n || throw(DimensionMismatch(
        "rp has length $(length(rp)) but linear solver has n=$(ls.n)"
    ))
    length(rd) == ls.m || throw(DimensionMismatch(
        "rd has length $(length(rd)) but linear solver has m=$(ls.m)"
    ))

    # Update diagonal scaling
    ls.θ .= θ
    # Update regularizers
    ls.rp .= rp
    ls.rd .= rd

    # Re-compute factorization
    # TODO: Keep S in memory, only change diagonal
    S = [
        spdiagm(0 => -(one(Tv) ./ ls.θ) .- rp)  ls.A';
        ls.A spdiagm(0 => rd)
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
- `ls::SparseIndefLinearSolver{Tv}`: Linear solver for the augmented system
- `A::AbstractMatrix{Tv}`: Constraint matrix
- `θ::Vector{Tv}`: Diagonal scaling vector
- `ξp::Vector{Tv}, ξd::Vector{Tv}`: Right-hand-side vectors
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::SparseIndefLinearSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
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
    # rp = A*dx + ls.rd .* dy - ξp
    # rd = - dx .* ( (one(Tv) ./ ls.θ) + ls.rp) + A' * dy - ξd
    # @info "\n|rp| = $(norm(rp, Inf))\n|rd| = $(norm(rd, Inf))\n"

    # ξ1 = [rd; rp]
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

"""
mutable struct SparsePosDefLinearSolver{Tv<:BlasReal, Ta<:SparseMatrixCSC{Tv, <:Integer}} <: PosDefLinearSolver{Tv, Ta}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::Ta
    θ::Vector{Tv}   # Diagonal scaling
    # Regularization
    rp::Vector{Tv}  # primal
    rd::Vector{Tv}  # dual

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

        return new{Tv, SparseMatrixCSC{Tv, Int}}(m, n, A, θ, zeros(Tv, n), ones(Tv, m), F)
    end

end

function update_linear_solver!(
    ls::SparsePosDefLinearSolver{Tv},
    θ::AbstractVector{Tv},
    rp::AbstractVector{Tv}=zeros(Tv, ls.n),
    rd::AbstractVector{Tv}=zeros(Tv, ls.m)
) where{Tv<:BlasReal}
    
    # Sanity checks
    length(θ) == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(rp) == ls.n || throw(DimensionMismatch(
        "rp has length $(length(rp)) but linear solver has n=$(ls.n)"
    ))
    length(rd) == ls.m || throw(DimensionMismatch(
        "rd has length $(length(rd)) but linear solver has m=$(ls.m)"
    ))

    ls.θ .= θ
    ls.rp .= rp  # Primal regularization is disabled for normal equations
    ls.rd .= rd

    # Re-compute factorization
    # D = (Θ^{-1} + Rp)^{-1}
    D = Diagonal(
        one(Tv) ./ ( (one(Tv) ./ ls.θ) .+ ls.rp)
    )
    Rd = spdiagm(0 => ls.rd)
    S = ls.A * D * ls.A' + Rd

    # Update factorization
    cholesky!(ls.F, Symmetric(S), check=false)
    issuccess(ls.F) || throw(PosDefException(0))

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
- `ls::SparseIndefLinearSolver{Tv}`: Linear solver for the augmented system
- `A::AbstractMatrix{Tv}`: Constraint matrix
- `θ::Vector{Tv}`: Diagonal scaling vector
- `ξp::Vector{Tv}, ξd::Vector{Tv}`: Right-hand-side vectors
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::SparsePosDefLinearSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    ξ_ = ξp .+ ls.A * (ξd ./ ( (one(Tv) ./ ls.θ) .+ ls.rp))

    # Solve augmented system
    dy .= (ls.F \ ξ_)

    # Recover dx
    dx .= (ls.A' * dy - ξd) ./ ( (one(Tv) ./ ls.θ) .+ ls.rp)

    # TODO: Iterative refinement
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    return nothing
end