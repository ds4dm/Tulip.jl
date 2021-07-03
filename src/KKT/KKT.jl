module KKT

using LinearAlgebra
using SparseArrays

using LinearOperators

export AbstractKKTSystem, AbstractKKTBackend, AbstractKKTSolver
export KKTOptions

"""
    AbstractKKTSystem

Abstract type for KKT systems
"""
abstract type AbstractKKTSystem end

include("systems.jl")

"""
    AbstractKKTBackend

Abstract type for KKT backend, i.e., the actual linear solver.
"""
abstract type AbstractKKTBackend end

"""
    DefaultKKTBackend

Default setting for KKT backend.

Currently defaults to [`CholmodBackend`](@ref) for `Float64` arithmetic,
    and [`LDLFactBackend`](@ref) otherwise.
"""
struct DefaultKKTBackend <: AbstractKKTBackend end

"""
    AbstractKKTSolver{T}

Abstract container for solving KKT systems in arithmetic `T`.
"""
abstract type AbstractKKTSolver{T} end

"""
    KKTOptions{T}

KKT solver options.
"""
Base.@kwdef mutable struct KKTOptions{T}
    Backend::AbstractKKTBackend = DefaultKKTBackend()
    System::AbstractKKTSystem = DefaultKKTSystem()
end

"""
    setup(A, system, backend; kwargs...)

Instantiate a KKT solver object.
"""
function setup end

#
# Specialized implementations should extend the functions below
#

"""
    update!(kkt, θinv, regP, regD)

Update internal data and factorization/pre-conditioner.

After this call, `kkt` can be used to solve the augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
for given right-hand sides `ξd` and `ξp`.

# Arguments
* `kkt::AbstractKKTSolver{T}`: the KKT solver object
* `θinv::AbstractVector{T}`: ``θ⁻¹``
* `regP::AbstractVector{T}`: primal regularizations
* `regD::AbstractVector{T}`: dual regularizations
"""
function update! end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the symmetric quasi-definite augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
and over-write `dx`, `dy` with the result.

# Arguments
- `dx, dy`: Vectors of unknowns, modified in-place
- `kkt`: Linear solver for the augmented system
- `ξp, ξd`: Right-hand-side vectors
"""
function solve! end

"""
    arithmetic(kkt::AbstractKKTSolver)

Return the arithmetic used by the solver.
"""
arithmetic(::AbstractKKTSolver{T}) where T = T

"""
    backend(kkt)

Return the name of the solver's backend.
"""
backend(::AbstractKKTSolver) = "Unkown"

"""
    linear_system(kkt)

Return which system is solved by the kkt solver.
"""
linear_system(::AbstractKKTSolver) = "Unkown"

# Generic tests
include("Test/test.jl")

# Custom linear solvers
include("Dense/lapack.jl")
include("Cholmod/cholmod.jl")
include("LDLFactorizations/ldlfact.jl")
const TlpLDLFact = TlpLDLFactorizations
include("krylov.jl")

# Default backend and system choices
function setup(A, ::DefaultKKTSystem, ::DefaultKKTBackend)
    T = eltype(A)
    if T == Float64
        return setup(A, K2(), TlpCholmod.Backend())
    else
        return setup(A, K2(), TlpLDLFact.Backend())
    end
end

end  # module
