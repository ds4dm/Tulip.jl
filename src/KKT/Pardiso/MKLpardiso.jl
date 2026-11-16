module TlpMKLPardiso

using LinearAlgebra
using SparseArrays

import Pardiso

using ..KKT: AbstractKKTBackend, AbstractKKTSolver
using ..KKT: AbstractKKTSystem, K1, K2
import ..KKT: setup, update!, solve!, backend, linear_system

"""
    Backend

MKLPardisoSQD backend for solving linear systems.

See [`MKLPardisoSQD`](@ref) for further details.
"""
struct Backend <: AbstractKKTBackend end

"""
    MKLPardisoSQD{T,S<:AbstractKKTSystem}

MKLPardiso-based KKT solver.

# Supported arithmetics
* `Float64`

# Supported systems
* [`K2`](@ref) via ``LDLᵀ`` factorization

# Examples

* To solve the augmented system with MKLPardiso's ``LDL^{T}`` factorization:
```julia
set_parameter(tlp_model, "KKT_Backend", Tulip.KKT.MKLPardisoSQD.Backend())
set_parameter(tlp_model, "KKT_System", Tulip.KKT.K2())
```

# Remark

The MKLPardiso solver works only on Intel/AMD platforms.
"""
mutable struct MKLPardisoSQD{S} <: AbstractKKTSolver{Float64}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # Problem data
    A::SparseMatrixCSC{Float64, Int}
    θ::Vector{Float64}
    regP::Vector{Float64}  # primal regularization
    regD::Vector{Float64}  # dual regularization

    # Left-hand side matrix
    K::SparseMatrixCSC{Float64, Int}
    ξ::Vector{Float64}                             # RHS of KKT system

    # Linear solver
    ps::Pardiso.MKLPardisoSolver
end

function setup(A::SparseMatrixCSC{Float64, Int}, ::K2, ::Backend)

    m, n = size(A)
    θ = ones(n)
    regP = ones(Float64, n)
    regD = ones(Float64, m)
    ξ = zeros(Float64, m+n)

    # We store we lower-triangular of the matrix
    K = [
        spdiagm(0 => -θ)  spzeros(n, m);
        A spdiagm(0 => ones(m))
    ]

    ps = Pardiso.MKLPardisoSolver()

    # We use symmetric indefinite matrices
    Pardiso.pardisoinit(ps)
    Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

    # Set number of threads
    Pardiso.set_nprocs!(ps, 1)

    # Do the analysis
    Pardiso.set_phase!(ps, Pardiso.ANALYSIS)
    Pardiso.pardiso(ps, K, ones(m+n))

    return MKLPardisoSQD{K2}(m, n, A, θ, regP, regD, K, ξ, ps)

end

backend(::MKLPardisoSQD) = "MKLPardiso"
linear_system(::MKLPardisoSQD) = "Augmented system (K2)"

function update!(kkt::MKLPardisoSQD{K2}, θ, regP, regD)
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
    # K is stored as lower-triangular and only its diagonal changes.
    @inbounds for j in 1:kkt.n
        k = kkt.K.colptr[j]
        kkt.K.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.K.colptr[kkt.n+i]
        kkt.K.nzval[k] = regD[i]
    end

    # Compute numerical factorization
    Pardiso.set_phase!(kkt.ps, Pardiso.NUM_FACT)
    Pardiso.pardiso(kkt.ps, kkt.K, zeros(kkt.m + kkt.n))

    return nothing
end

function solve!(dx, dy, kkt::MKLPardisoSQD{K2}, ξp, ξd)
    m, n = kkt.m, kkt.n

    # Setup right-hand side
    kkt.ξ = [ξd; ξp]
#    @views copyto!(kkt.ξ[1:n], ξd)
#    @views copyto!(kkt.ξ[(n+1):end], ξp)

    # Solve augmented system
    d = zeros(m + n)
    Pardiso.set_phase!(kkt.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    Pardiso.pardiso(kkt.ps, d, kkt.K, kkt.ξ)

    # Recover dx, dy
#    @views copyto!(dx, kkt.ξ[1:n])
#    @views copyto!(dy, kkt.ξ[(n+1):end])
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # TODO: iterative refinement
    return nothing
end

end  # module
