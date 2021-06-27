using LinearAlgebra:BlasReal

"""
    DenseBackend

Dense linear algebra backend for solving linear systems.
"""
struct DenseBackend <: AbstractKKTBackend end

"""
    DenseSolver{T,S}

Dense linear algebra-based KKT solver.

Instantiated through [`DenseBackend`](@ref).

# Supported arithmetics

All arithmetics are supported.
BLAS/LAPACK routines are used automatically with `Float32` and `Float64` arithmetic.

# Supported systems
* [`K1`](@ref) via Cholesky factorization
"""
mutable struct DenseSolver{T} <: AbstractKKTSolver{T}
    # Problem data
    m::Int
    n::Int
    A::Matrix{T}

    # Workspace
    _A::Matrix{T}    # Place-holder for scaled copy of A
    θ::Vector{T}     # Diagonal scaling
    regP::Vector{T}  # Primal regularization
    regD::Vector{T}  # Dual regularization
    K::Matrix{T}     # KKT matrix
    ξ::Vector{T}     # RHS of KKT system
end

backend(::DenseSolver) = "Julia.LinearAlgebra"
backend(::DenseSolver{<:BlasReal}) = "LAPACK $(LinearAlgebra.BLAS.vendor())"
linear_system(::DenseSolver) = "Normal equations (K1)"

function setup(A::Matrix{T}, ::K1, ::DenseBackend) where{T}
    m, n = size(A)

    _A = Matrix{T}(undef, m, n)
    θ = ones(T, n)
    regP = ones(T, n)
    regD = ones(T, m)
    K = Matrix{T}(undef, m, m)
    ξ = zeros(T, m)

    return DenseSolver{T}(m, n, A, _A, θ, regP, regD, K, ξ)
end

function update!(kkt::DenseSolver, θ, regP, regD)
    m, n = kkt.m, kkt.n

    # Sanity checks
    length(θ)  == n || throw(DimensionMismatch(
        "length(θ)=$(length(θ)) but KKT solver has n=$n."
    ))
    length(regP) == n || throw(DimensionMismatch(
        "length(regP)=$(length(regP)) but KKT solver has n=$n"
    ))
    length(regD) == m || throw(DimensionMismatch(
        "length(regD)=$(length(regD)) but KKT solver has m=$m"
    ))

    copyto!(kkt.θ, θ)
    copyto!(kkt.regP, regP)
    copyto!(kkt.regD, regD)

    # Re-compute normal equations matrix
    # There's no function that does S = A*D*A', so we cache a copy of A
    copyto!(kkt._A, kkt.A)
    D = sqrt(inv(Diagonal(kkt.θ .+ kkt.regP)))
    rmul!(kkt._A, D)  # B = A * √D
    mul!(kkt.K, kkt._A, transpose(kkt._A), true, false)  # Now K = A*D*A'
    # Finally, add dual regularizations to the diagonal
    @inbounds for i in 1:kkt.m
        kkt.K[i, i] += kkt.regD[i]
    end

    # In-place Cholesky factorization
    cholesky!(Symmetric(kkt.K))

    return nothing
end

function solve!(dx, dy, kkt::DenseSolver{T}, ξp, ξd) where{T}
    m, n = kkt.m, kkt.n

    # Set-up right-hand side
    D = inv(Diagonal(kkt.θ .+ kkt.regP))
    copyto!(dy, ξp)
    mul!(dy, kkt.A, D * ξd, true, true)

    # Solve normal equations
    ldiv!(UpperTriangular(kkt.K)', dy)
    ldiv!(UpperTriangular(kkt.K) , dy)

    # Recover dx
    copyto!(dx, ξd)
    mul!(dx, kkt.A', dy, one(T), -one(T))
    lmul!(D, dx)

    # TODO: Iterative refinement
    return nothing
end
