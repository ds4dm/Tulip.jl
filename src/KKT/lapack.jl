using LinearAlgebra.LAPACK

const BlasReal = LinearAlgebra.BlasReal

"""
    DenseSPD{T}

KKT solver for dense linear systems.

Uses a Cholesky factorization of the normal equations system. Not recommended
    for systems of size larger than a few thousands.

```julia
model = Tulip.Model{Float64}()
model.params.KKTOptions = Tulip.KKT.SolverOptions(Tulip.KKT.DenseSPD)
```

!!! info
    Dispatches to BLAS/LAPACK in `Float32` and `Float64` arithmetic.
"""
mutable struct DenseSPD{T} <: AbstractKKTSolver{T}
    m::Int  # Number of rows in A
    n::Int  # Number of columns in A

    # Problem data
    A::Matrix{T}
    θ::Vector{T}

    # Regularizations
    regP::Vector{T}  # primal
    regD::Vector{T}  # dual

    _A::Matrix{T}  # place-holder to avoid allocating memory afterwards

    # Factorization
    S::Matrix{T}  # Normal equations matrix, that also holds the factorization

    # Constructor
    function DenseSPD(A::AbstractMatrix{T}) where{T}
        m, n = size(A)
        θ = ones(T, n)

        # We just need to allocate memory here,
        # so no factorization yet
        S = zeros(T, m, m)

        return new{T}(m, n, A, θ, zeros(T, n), zeros(T, m), Matrix(A), S)
    end
end

backend(::DenseSPD) = "LAPACK"
linear_system(::DenseSPD) = "Normal equations"

"""
    update!(kkt::DenseSPD{T}, θ, regP, regD)

Compute normal equations system matrix and update Cholesky factorization.

Uses Julia's generic linear algebra.
"""
function update!(
    kkt::DenseSPD{T},
    θ::AbstractVector{T},
    regP::AbstractVector{T},
    regD::AbstractVector{T}
) where{T}
    # Sanity checks
    length(θ) == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    kkt.θ .= θ
    kkt.regP .= regP
    kkt.regD .= regD

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(kkt._A, kkt.A)
    D = Diagonal(one(T) ./ (kkt.θ .+ kkt.regP))
    rmul!(kkt._A, D)  # B = A * D
    mul!(kkt.S, kkt._A, transpose(kkt.A))  # Now S = A*D*A'
    # TODO: do not re-compute S if only dual regularization changes
    @inbounds for i in 1:kkt.m
        kkt.S[i, i] += kkt.regD[i]
    end
    
    # Cholesky factorization
    cholesky!(Symmetric(kkt.S))

    return nothing
end

"""
    update!(kkt::DenseSPD{<:BlasReal}, θ, regP, regD)

Compute normal equations system matrix and update Cholesky factorization.

Uses BLAS and LAPACK routines.
"""
function update!(
    kkt::DenseSPD{T},
    θ::AbstractVector{T},
    regP::AbstractVector{T},
    regD::AbstractVector{T}
) where{T<:BlasReal}
    # Sanity checks
    length(θ) == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    kkt.θ .= θ
    kkt.regP .= regP
    kkt.regD .= regD

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(kkt._A, kkt.A)
    D = Diagonal(sqrt.(one(T) ./ (kkt.θ .+ kkt.regP)))
    rmul!(kkt._A, D)  # B = A * √D
    BLAS.syrk!('U', 'N', one(T), kkt._A, zero(T), kkt.S)

    # Regularization
    # TODO: only update diagonal of S if only dual regularization changes
    @inbounds for i in 1:kkt.m
        kkt.S[i, i] += kkt.regD[i]
    end
    
    # Cholesky factorization
    # TODO: regularize if needed
    LAPACK.potrf!('U', kkt.S)
    
    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.

Uses two generic triangular solves for solving the normal equations system.
"""
function solve!(
    dx::Vector{T}, dy::Vector{T},
    kkt::DenseSPD{T},
    ξp::Vector{T}, ξd::Vector{T}
) where{T}
    m, n = kkt.m, kkt.n
    
    # Set-up right-hand side
    dy .= ξp .+ kkt.A * (ξd ./ (kkt.θ .+ kkt.regP))

    # Solve normal equations
    ldiv!(UpperTriangular(kkt.S)', dy)
    ldiv!(UpperTriangular(kkt.S) , dy)

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (kkt.A' * dy - ξd) ./ (kkt.θ .+ kkt.regP)

    # TODO: Iterative refinement
    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.

Uses one LAPACK call for solving the normal equations system.
"""
function solve!(
    dx::Vector{T}, dy::Vector{T},
    kkt::DenseSPD{T},
    ξp::Vector{T}, ξd::Vector{T}
) where{T<:BlasReal}
    m, n = kkt.m, kkt.n
    
    # Set-up right-hand side
    dy .= ξp .+ kkt.A * (ξd ./ (kkt.θ .+ kkt.regP))

    # Solve augmented system
    LAPACK.potrs!('U', kkt.S, dy)  # potrs! over-writes the right-hand side

    # Recover dx
    # TODO: use more efficient mul! syntax
    dx .= (kkt.A' * dy - ξd) ./ (kkt.θ .+ kkt.regP)

    # TODO: Iterative refinement
    return nothing
end