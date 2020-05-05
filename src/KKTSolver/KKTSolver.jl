module KKT

using LinearAlgebra

export AbstractKKTSolver

"""
    AbstractKKTSolver{Tv}

Abstract container for solving an augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
where `ξd` and `ξp` are given right-hand side.
"""
abstract type AbstractKKTSolver{Tv} end

"""
    SolverOptions

Used to pass options and instantiate KKT solvers.
"""
struct SolverOptions
    Ts::Type{<:AbstractKKTSolver}
    options::Base.Iterators.Pairs

    SolverOptions(::Type{Ts}; kwargs...) where{Ts<:AbstractKKTSolver} = new(Ts, kwargs)
end

"""
    setup(T::Type{<:AbstractKKTSolver}, args...; kwargs...)

Instantiate a KKT solver object.
"""
function setup(Ts::Type{<:AbstractKKTSolver}, args...; kwargs...)
    return Ts(args...; kwargs...)
end

# 
# Specialized implementations should extend the functions below
# 

"""
    update!(ls, θinv, regP, regD)

Update internal data and factorization/pre-conditioner.

After this call, `ls` can be used to solve the augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
for given right-hand sides `ξd` and `ξp`.

# Arguments
* `ls::AbstractKKTSolver{Tv}`: the KKT solver object
* `θinv::AbstractVector{Tv}`: ``θ⁻¹``
* `regP::AbstractVector{Tv}`: primal regularizations
* `regD::AbstractVector{Tv}`: dual regularizations
"""
function update! end


"""
    solve!(dx, dy, ls, ξp, ξd)

Solve the symmetric quasi-definite augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
and over-write `dx`, `dy` with the result.

# Arguments
- `dx, dy`: Vectors of unknowns, modified in-place
- `ls`: Linear solver for the augmented system
- `ξp, ξd`: Right-hand-side vectors
"""
function solve! end

"""
    arithmetic(kkt::AbstractKKTSolver)

Return the arithmetic used by the solver.
"""
arithmetic(kkt::AbstractKKTSolver{Tv}) where Tv = Tv

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
include("test.jl")

# Custom linear solvers
include("lapack.jl")
include("cholmod.jl")
include("ldlfact.jl")

"""
    default_options(Tv)

Use CHOLMOD for `Float64` and LDLFactorizations otherwise.
"""
function default_options(::Type{Tv}) where{Tv}
    if Tv == Float64
        return SolverOptions(CholmodSolver; normal_equations=false)
    else
        return SolverOptions(LDLFact_SymQuasDef)
    end
end

end  # module