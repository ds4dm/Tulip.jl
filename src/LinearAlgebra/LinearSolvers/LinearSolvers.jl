"""
    AbstractLinearSolver{Tv}

Abstract container for linear solver used in solving the augmented system.
"""
abstract type AbstractLinearSolver{Tv<:Real} end

#= TODO: Use traits instead, e.g.
    * Direct vs indirect method
    * Numerical precision (Tv)
    * Augmented system vs Normal equations systems
=#
"""
    IndefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the indefinite augmented system.
"""
abstract type IndefLinearSolver{Tv<:Real} <: AbstractLinearSolver{Tv} end

"""
    PosDefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the PSD normal equations system.
"""
abstract type PosDefLinearSolver{Tv<:Real} <: AbstractLinearSolver{Tv} end


# 
# Specialized implementations should extend the two functions below.
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


# 
include("dense.jl")
include("sparse.jl")
include("LDLF.jl")

# TODO: use parameter to choose between Indef/PosDef system
AbstractLinearSolver(A::Matrix{Tv}) where{Tv<:Real} = DenseLinearSolver(A)
AbstractLinearSolver(
    A::SparseMatrixCSC{Tv, Int64}
) where{Tv<:BlasReal} = SparseIndefLinearSolver(A)
AbstractLinearSolver(
    A::SparseMatrixCSC{Tv, Int64}
) where{Tv<:Real} = LDLFLinearSolver(A)