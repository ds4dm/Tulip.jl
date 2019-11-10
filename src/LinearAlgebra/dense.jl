using LinearAlgebra.LAPACK

function construct_matrix(
    ::Type{Matrix}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{Tv}
) where{Tv<:Real}
    A = zeros(Tv, m, n)
    # TODO: may be more efficient to first sort indices so that
    # A is accessed in column-major order.
    for(i, j, v) in zip(aI, aJ, aV)
        A[i, j] = v
    end
    return A
end

"""
    DenseLinearSolver{Tv}

Linear solver for dense matrices.

The augmented system is automatically reduced to the normal equations system.
BLAS/LAPACK functions are used whenever applicable.
"""
mutable struct DenseLinearSolver{Tv<:Real} <: TLPLinearSolver{Tv}
    m::Int  # Number of rows in A
    n::Int  # Number of columns in A

    # Problem data
    A::Matrix{Tv}
    θ::Vector{Tv}

    # Regularizations
    rp::Vector{Tv}  # primal
    rd::Vector{Tv}  # dual

    B::Matrix{Tv}  # place-holder to avoid allocating memory afterwards

    # Factorization
    S::Matrix{Tv}  # Normal equations matrix, that also holds the factorization

    # Constructor
    function DenseLinearSolver(A::Matrix{Tv}) where{Tv<:Real}
        m, n = size(A)
        θ = ones(Tv, n)

        # We just need to allocate memory here,
        # so no factorization yet
        S = zeros(Tv, m, m)

        return new{Tv}(m, n, A, θ, zeros(Tv, n), zeros(Tv, m), copy(A), S)
    end
end

TLPLinearSolver(A::Matrix{Tv}) where{Tv<:Real} = DenseLinearSolver(A)

# generic
function update_linear_solver(
    ls::DenseLinearSolver{Tv},
    d::AbstractVector{Tv},
    rp::AbstractVector{Tv}=zeros(Tv, ls.n),
    rd::AbstractVector{Tv}=zeros(Tv, ls.m)
) where{Tv<:Real}
    # Sanity checks
    length(d) == ls.n || throw(DimensionMismatch(
        "d has length $(length(d)) but linear solver is for n=$(ls.n)."
    ))
    length(rp) == ls.n || throw(DimensionMismatch(
        "rp has length $(length(rp)) but linear solver has n=$(ls.n)"
    ))
    length(rd) == ls.m || throw(DimensionMismatch(
        "rd has length $(length(rd)) but linear solver has m=$(ls.m)"
    ))

    ls.θ .= d
    ls.rp .= rp
    ls.rd .= rd

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(ls.B, ls.A)
    D = Diagonal(
        one(Tv) ./ ( (one(Tv) ./ ls.θ) .+ ls.rp)
    )
    rmul!(ls.B, D)  # B = A * D
    mul!(ls.S, ls.B, transpose(ls.A))  # Now S = A*D*A'
    # TODO: do not re-compute S if only dual regularization changes
    @inbounds for i in 1:ls.m
        ls.S[i, i] += ls.rd[i]
    end
    
    # Cholesky factorization
    cholesky!(Symmetric(ls.S))

    return nothing
end

# Use BLAS/LAPACK if available
function update_linear_solver(
    ls::DenseLinearSolver{Tv},
    d::AbstractVector{Tv},
    rp::AbstractVector{Tv}=zeros(Tv, ls.n),
    rd::AbstractVector{Tv}=zeros(Tv, ls.m)
) where{Tv<:BlasReal}
    # Sanity checks
    length(d) == ls.n || throw(DimensionMismatch(
        "d has length $(length(d)) but linear solver is for n=$(ls.n)."
    ))
    length(rp) == ls.n || throw(DimensionMismatch(
        "rp has length $(length(rp)) but linear solver has n=$(ls.n)"
    ))
    length(rd) == ls.m || throw(DimensionMismatch(
        "rd has length $(length(rd)) but linear solver has m=$(ls.m)"
    ))

    ls.θ .= d
    ls.rp .= rp
    ls.rd .= rd

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(ls.B, ls.A)
    D = Diagonal(sqrt.(
        one(Tv) ./ ( (one(Tv) ./ ls.θ) .+ ls.rp)
    ))
    rmul!(ls.B, D)  # B = A * √D
    BLAS.syrk!('U', 'N', one(Tv), ls.B, zero(Tv), ls.S)

    # Regularization
    # TODO: only update diagonal of S if only dual regularization changes
    @inbounds for i in 1:ls.m
        ls.S[i, i] += ls.rd[i]
    end
    
    # Cholesky factorization
    # TODO: regularize if needed
    LAPACK.potrf!('U', ls.S)
    
    return nothing
end


"""
    solve_augmented_system!(dx, dy, A, ls, θ, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
```
    -Θ*dx  + A'dy = ξd
     A*dx         = ξp
```

For dense `A` we first compute the normal equations system.

`A` and `θ` are given to perform iterative refinement if needed.

# Arguments
- `dx::Vector{Tv}, dy::Vector{Tv}`: Vectors of unknowns, modified in-place
- `ls::DenseLinearSolver{Tv}`: Linear solver for the augmented system
- `A::Matrix{Tv}`: Constraint matrix
- `θ::Vector{Tv}`: Diagonal scaling vector
- `ξp::Vector{Tv}, ξd::Vector{Tv}`: Right-hand-side vectors
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::DenseLinearSolver{Tv}, A::Matrix{Tv}, θ::Vector{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = size(A)
    
    # Set-up right-hand side
    dy .= ξp .+ ls.A * (ξd ./ ( (one(Tv) ./ θ) .+ ls.rp))

    # Solve augmented system
    # 
    ldiv!(UpperTriangular(ls.S)', dy)
    ldiv!(UpperTriangular(ls.S) , dy)

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (ls.A' * dy - ξd) ./ ( (one(Tv) ./ θ) .+ ls.rp)

    # TODO: Iterative refinement
    return nothing
end

function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::DenseLinearSolver{Tv}, A::Matrix{Tv}, θ::Vector{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:BlasReal}
    m, n = size(A)
    
    # Set-up right-hand side
    dy .= ξp .+ ls.A * (ξd ./ ( (one(Tv) ./ θ) .+ ls.rp))

    # Solve augmented system
    LAPACK.potrs!('U', ls.S, dy)  # potrs! over-writes the right-hand side

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (ls.A' * dy - ξd) ./ ( (one(Tv) ./ θ) .+ ls.rp)

    # TODO: Iterative refinement
    return nothing
end