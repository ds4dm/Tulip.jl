module TlpCholmod

using LinearAlgebra
using SparseArrays
using SuiteSparse.CHOLMOD

using ..KKT: AbstractKKTBackend, AbstractKKTSolver
using ..KKT: AbstractKKTSystem, K1, K2
import ..KKT: setup, update!, solve!, backend, linear_system

"""
    Backend

CHOLMOD backend for solving linear systems.

See [`CholmodSolver`](@ref) for further details.
"""
struct Backend <: AbstractKKTBackend end

"""
    CholmodSolver{T,S<:AbstractKKTSystem}

CHOLMOD-based KKT solver.

# Supported arithmetics
* `Float64`

# Supported systems
* [`K2`](@ref) via ``LDLᵀ`` factorization
* [`K1`](@ref) via Cholesky (``LLᵀ``) factorization

# Examples

* To solve the augmented system with CHOLMOD's ``LDL^{T}`` factorization:
```julia
set_parameter(tlp_model, "KKT_Backend", Tulip.KKT.TlpCholmod.Backend())
set_parameter(tlp_model, "KKT_System", Tulip.KKT.K2())
```

* To solve the normal equations system with CHOLMOD's Cholesky factorization:
```julia
set_parameter(tlp_model, "KKT_Backend", Tulip.KKT.TlpCholmod.Backend())
set_parameter(tlp_model, "KKT_System", Tulip.KKT.K1())
```
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
