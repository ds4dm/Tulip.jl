using SuiteSparse.CHOLMOD

"""
    CholmodBackend

CHOLMOD backend for solving linear systems.
"""
struct CholmodBackend <: AbstractKKTBackend end

"""
    CholmodSolver{T,S}

CHOLMOD-based KKT solver for system `S` and arithmetic `T`.

Instantiated by using [`CholmodBackend`](@ref).

# Supported arithmetics
* `Float64`

# Supported systems
* [`K2`](@ref) via ``L×D×Lᵀ`` factorization
* [`K1`](@ref) via Cholesky (``L×Lᵀ``) factorization
"""
mutable struct CholmodSolver{T,S} <: AbstractKKTSolver{T}
    # Problem data
    m::Int
    n::Int
    A::SparseMatrixCSC{T,Int}

    # Workspace
    # TODO: store K as CHOLMOD.Sparse instead of SparseMatrixCSC
    θ::Vector{T}               # Diagonal scaling
    regP::Vector{T}            # Primal regularization
    regD::Vector{T}            # Dual regularization
    K::SparseMatrixCSC{T,Int}  # KKT matrix
    F::CHOLMOD.Factor{T}       # Factorization
    ξ::Vector{T}               # RHS of KKT system
end

const CholmodSQD = CholmodSolver{Float64,K2}
const CholmodSPD = CholmodSolver{Float64,K1}

backend(::CholmodSolver) = "CHOLMOD"
linear_system(::CholmodSQD) = "Augmented system (K2)"
linear_system(::CholmodSPD) = "Normal equations (K1)"

# Convert to sparse matrix if other type is used
setup(A, system, backend::CholmodBackend) = setup(convert(SparseMatrixCSC, A), system, backend)

function setup(A::SparseMatrixCSC{Float64,Int}, ::K2, ::CholmodBackend)
    m, n = size(A)

    θ = ones(Float64, n)
    regP = ones(Float64, n)
    regD = ones(Float64, m)
    ξ = zeros(Float64, m+n)

    K = [
        spdiagm(0 => -θ)  A';
        spzeros(Float64, m, n) spdiagm(0 => ones(m))
    ]

    # TODO: Symbolic factorization only
    F = ldlt(Symmetric(K))

    return CholmodSolver{Float64,K2}(m, n, A, θ, regP, regD, K, F, ξ)
end

function setup(A::SparseMatrixCSC{Float64}, ::K1, ::CholmodBackend)
    m, n = size(A)

    θ = ones(Float64, n)
    regP = ones(Float64, n)
    regD = ones(Float64, m)
    ξ = zeros(Float64, m)

    # TODO: analyze + in-place A*D*A' product
    K = sparse(A * A') + spdiagm(0 => regD)

    # TODO: PSD-ness checks
    F = cholesky(Symmetric(K))

    return CholmodSolver{Float64,K1}(m, n, A, θ, regP, regD, K, F, ξ)
end

function update!(kkt::CholmodSQD, θ, regP, regD)
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

    # Update KKT matrix
    # K is stored as upper-triangular, and only its diagonal is changed
    @inbounds for j in 1:kkt.n
        k = kkt.K.colptr[1+j] - 1
        kkt.K.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.K.colptr[1+kkt.n+i] - 1
        kkt.K.nzval[k] = regD[i]
    end

    ldlt!(kkt.F, Symmetric(kkt.K))
    return nothing
end

function update!(kkt::CholmodSPD, θ, regP, regD)
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

    # Form normal equations matrix
    # TODO: use in-place update of S
    D = inv(Diagonal(kkt.θ .+ kkt.regP))
    kkt.K = (kkt.A * D * kkt.A') + spdiagm(0 => kkt.regD)

    # Update factorization
    cholesky!(kkt.F, Symmetric(kkt.K), check=false)
    issuccess(kkt.F) || throw(PosDefException(0))

    return nothing
end

function solve!(dx, dy, kkt::CholmodSQD, ξp, ξd)
    m, n = kkt.m, kkt.n

    # Setup right-hand side
    @views copyto!(kkt.ξ[1:n], ξd)
    @views copyto!(kkt.ξ[(n+1):end], ξp)

    # Solve augmented system
    # CHOLMOD doesn't have in-place solve, so this line will allocate
    δ = kkt.F \ kkt.ξ

    # Recover dx, dy
    @views copyto!(dx, δ[1:n])
    @views copyto!(dy, δ[(n+1):end])

    # TODO: iterative refinement
    return nothing
end

function solve!(dx, dy, kkt::CholmodSPD, ξp, ξd)
    m, n = kkt.m, kkt.n

    D = inv(Diagonal(kkt.θ .+ kkt.regP))
    copyto!(kkt.ξ, ξp)
    mul!(kkt.ξ, kkt.A, D * ξd, true, true)

    # Solve normal equations
    # CHOLMOD doesn't have in-place solve, so this line will allocate
    dy .= kkt.F \ kkt.ξ

    # Recover dx
    copyto!(dx, ξd)
    mul!(dx, kkt.A', dy, 1.0, -1.0)
    lmul!(D, dx)

    # TODO: iterative refinement
    return nothing
end
