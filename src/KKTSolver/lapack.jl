using LinearAlgebra.LAPACK

"""
    Lapack <: LSBackend

Use LAPACK backend.

Options available:
* Numerical precision: `Tv<:Real`
    * `Float32` and `Float64` will use LAPACK.
    * Other numerical types will use Julia's generic cholesky factorization.
* [`NormalEquations`](@ref) only, with Cholesky factorization
"""
struct Lapack <: LSBackend end

"""
    KKTSolver_Dense{Tv}

Linear solver for the 2x2 augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [xi_d]
    [   A          Rd] [dy]   [xi_p]
```
with ``A`` dense.

The augmented system is automatically reduced to the normal equations system.
BLAS/LAPACK functions are used whenever applicable.
"""
mutable struct KKTSolver_Dense{Tv<:Real} <: AbstractKKTSolver{Tv}
    m::Int  # Number of rows in A
    n::Int  # Number of columns in A

    # Problem data
    A::Matrix{Tv}
    θ::Vector{Tv}

    # Regularizations
    regP::Vector{Tv}  # primal
    regD::Vector{Tv}  # dual

    B::Matrix{Tv}  # place-holder to avoid allocating memory afterwards

    # Factorization
    S::Matrix{Tv}  # Normal equations matrix, that also holds the factorization

    # Constructor
    function KKTSolver_Dense(A::Matrix{Tv}) where{Tv<:Real}
        m, n = size(A)
        θ = ones(Tv, n)

        # We just need to allocate memory here,
        # so no factorization yet
        S = zeros(Tv, m, m)

        return new{Tv}(m, n, A, θ, zeros(Tv, n), zeros(Tv, m), copy(A), S)
    end
end

AbstractKKTSolver(
    ::Lapack,
    ::NormalEquations,
    A::AbstractMatrix{Tv}
) where{Tv<:Real} = KKTSolver_Dense(Matrix(A))

backend(::KKTSolver_Dense) = "LAPACK"
linear_system(::KKTSolver_Dense) = "Normal equations"

"""
    update_linear_solver!(ls::KKTSolver_Dense{<:Real}, θ, regP, regD)

Compute normal equations system matrix and update Cholesky factorization.

Uses Julia's generic linear algebra.
"""
function update_linear_solver!(
    ls::KKTSolver_Dense{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:Real}
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
    ls.regP .= regP
    ls.regD .= regD

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(ls.B, ls.A)
    D = Diagonal(one(Tv) ./ (ls.θ .+ ls.regP))
    rmul!(ls.B, D)  # B = A * D
    mul!(ls.S, ls.B, transpose(ls.A))  # Now S = A*D*A'
    # TODO: do not re-compute S if only dual regularization changes
    @inbounds for i in 1:ls.m
        ls.S[i, i] += ls.regD[i]
    end
    
    # Cholesky factorization
    cholesky!(Symmetric(ls.S))

    return nothing
end

"""
    update_linear_solver!(ls::KKTSolver_Dense{<:BlasReal}, θ, regP, regD)

Compute normal equations system matrix and update Cholesky factorization.

Uses BLAS and LAPACK routines.
"""
function update_linear_solver!(
    ls::KKTSolver_Dense{Tv},
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
    ls.regP .= regP
    ls.regD .= regD

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(ls.B, ls.A)
    D = Diagonal(sqrt.(one(Tv) ./ (ls.θ .+ ls.regP)))
    rmul!(ls.B, D)  # B = A * √D
    BLAS.syrk!('U', 'N', one(Tv), ls.B, zero(Tv), ls.S)

    # Regularization
    # TODO: only update diagonal of S if only dual regularization changes
    @inbounds for i in 1:ls.m
        ls.S[i, i] += ls.regD[i]
    end
    
    # Cholesky factorization
    # TODO: regularize if needed
    LAPACK.potrf!('U', ls.S)
    
    return nothing
end


"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.

Uses two generic triangular solves for solving the normal equations system.
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::KKTSolver_Dense{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    dy .= ξp .+ ls.A * (ξd ./ (ls.θ .+ ls.regP))

    # Solve normal equations
    ldiv!(UpperTriangular(ls.S)', dy)
    ldiv!(UpperTriangular(ls.S) , dy)

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (ls.A' * dy - ξd) ./ (ls.θ .+ ls.regP)

    # TODO: Iterative refinement
    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.

Uses one LAPACK call for solving the normal equations system.
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::KKTSolver_Dense{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:BlasReal}
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    dy .= ξp .+ ls.A * (ξd ./ (ls.θ .+ ls.regP))

    # Solve augmented system
    LAPACK.potrs!('U', ls.S, dy)  # potrs! over-writes the right-hand side

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (ls.A' * dy - ξd) ./ (ls.θ .+ ls.regP)

    # TODO: Iterative refinement
    return nothing
end
