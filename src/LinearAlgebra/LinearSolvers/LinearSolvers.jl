"""
    AbstractLinearSolver{Tv}

Abstract container for solving an augmented system
```
    [-(Θ^{-1} + Rp)   A'] [dx] = [ξd]
    [   A             Rd] [dy] = [ξp]
```
where `ξd` and `ξp` are given right-hand side.
"""
abstract type AbstractLinearSolver{Tv<:Real} end


"""
    LinearSystem

Indicates which linear system is solved
* `DefaultSystem`: default 
* `AugmentedSystem`: solve the symmetric quasi-definite augmented system
* `NormalEquations`: solve the symmetric positive-definite normal equations 
"""
abstract type LinearSystem end

struct DefaultSystem <: LinearSystem end
struct AugmentedSystem <: LinearSystem end
struct NormalEquations <: LinearSystem end

"""
    LSBackend

Backend used for solving linear systems.
"""
abstract type LSBackend end

struct DefaultBackend <: LSBackend end

# 
# Specialized implementations should extend the functions below
# 

"""
    update_linear_solver!(ls, θ, regP, regD)

Update internal data, and re-compute factorization/pre-conditioner.

After this call, `ls` can be used to solve the augmented system
```
    [-(Θ^{-1} + Rp)   A'] [dx] = [ξd]
    [   A             Rd] [dy] = [ξp]
```
for given right-hand sides `ξd` and `ξp`.
"""
function update_linear_solver! end


"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the symmetric quasi-definite augmented system
```
    [-(Θ^{-1} + Rp)   A'] [dx] = [ξd]
    [   A             Rd] [dy] = [ξp]
```
and over-write `dx`, `dy` with the result.

# Arguments
- `dx, dy`: Vectors of unknowns, modified in-place
- `ls`: Linear solver for the augmented system
- `ξp, ξd`: Right-hand-side vectors
"""
function solve_augmented_system! end


include("test.jl")

# 
include("lapack.jl")
include("cholmod.jl")
include("ldlfact.jl")

# Default settings
AbstractLinearSolver(
    ::DefaultBackend,
    ::DefaultSystem,
    A::Matrix{Tv}
) where{Tv<:Real} = AbstractLinearSolver(Lapack(), NormalEquations(), A)

AbstractLinearSolver(
    ::DefaultBackend,
    ::DefaultSystem,
    A::AbstractMatrix{Float64}
) = AbstractLinearSolver(Cholmod(), AugmentedSystem(), A)

AbstractLinearSolver(
    ::DefaultBackend,
    ::DefaultSystem,
    A::AbstractMatrix{Tv}
) where{Tv<:Real} = AbstractLinearSolver(LDLFact(), AugmentedSystem(), A)