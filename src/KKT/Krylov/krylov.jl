module TlpKrylov

using LinearAlgebra

using Krylov
using LinearOperators
const LO = LinearOperators

using ..KKT: AbstractKKTBackend, AbstractKKTSolver
using ..KKT: AbstractKKTSystem, K1, K2
import ..KKT: arithmetic, backend, linear_system
import ..KKT: setup, update!, solve!

include("defs.jl")

"""
    Backend{KS<:Krylov.KrylovWorkspace,V<:AbstractVector}

[Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl)-based backend for solving linear systems.

The type is parametrized by:
* `KS<:Krylov.KrylovWorkspace`: workspace type for the Krylov method.
    Also defines the Krylov method to be used.
* `V<:AbstractVector`: the vector storage type used within the Krylov method.
    This should be set to `Vector{T}` (for arithmetic `T`) unless, e.g., one uses a GPU.

See the [Krylov.jl documentation](https://juliasmoothoptimizers.github.io/Krylov.jl/dev/inplace/) for further details.

# Example usage

All the following examples assume everything runs on a CPU in `Float64` arithmetic.
* To use the conjugate gradient:
```julia
backend = KKT.TlpKrylov.Backend(Krylov.CgWorkspace, Vector{Float64})
```
* To use MINRES:
```julia
backend = KKT.TlpKrylov.Backend(Krylov.MinresWorkspace, Vector{Float64})
```
"""
struct Backend{KS,V} <: AbstractKKTBackend
    krylov_solver::Type{KS}
    vector_storage::Type{V}
end

"""
    AbstractKrylovSolver{T}

Abstract type for Kyrlov-based linear solvers.
"""
abstract type AbstractKrylovSolver{T} <: AbstractKKTSolver{T} end

include("spd.jl")
include("sid.jl")
include("sqd.jl")

end  # module
