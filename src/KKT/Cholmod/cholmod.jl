module TlpCholmod

using LinearAlgebra
using SparseArrays
using SuiteSparse.CHOLMOD

using ..KKT: AbstractKKTBackend, AbstractKKTSolver
using ..KKT: AbstractKKTSystem, K1, K2
import ..KKT: setup, update!, solve!, backend, linear_system

"""
    Cholmod.Backend

CHOLMOD backend for solving linear systems.
"""
struct Backend <: AbstractKKTBackend end

"""
    Cholmod.KKTSolver{T,S}

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

backend(::CholmodSolver) = "CHOLMOD"

# Convert to sparse matrix if other type is used
setup(A, system, backend::Backend) = setup(convert(SparseMatrixCSC, A), system, backend)

include("spd.jl")  # Normal equations
include("sqd.jl")  # Augmented system

end  # module
