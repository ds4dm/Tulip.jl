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
    Backend{KS}

[Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl)-based backend for solving linear systems.

```julia
backend = KKT.TlpKrylov.Backend(CgSolver, Vector{Float64})
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

end  # module
