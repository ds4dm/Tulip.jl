"""
    AbstractLinearSolver{Tv, Ta}

Abstract container for linear solver used in solving the augmented system.
"""
abstract type AbstractLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} end

#= TODO: Use traits instead, e.g.
    * Direct vs indirect method
    * Numerical precision (Tv)
    * Augmented system vs Normal equations systems
=#
"""
    IndefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the indefinite augmented system.
"""
abstract type IndefLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} <: AbstractLinearSolver{Tv, Ta} end

"""
    PosDefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the PSD normal equations system.
"""
abstract type PosDefLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} <: AbstractLinearSolver{Tv, Ta} end


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
function update_linear_solver!(ls::AbstractLinearSolver, θ, regP, regD) end


"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the symmetric quasi-definite augmented system
```
    [-(Θ^{-1} + Rp)   A'] [dx] = [ξd]
    [   A             Rd] [dy] = [ξp]
```
and over-write `dx`, `dy` with the result.
"""
function solve_augmented_system!(
    dx, dy,
    ls::AbstractLinearSolver,
    ξp, ξd
) end


# 
include("dense.jl")
include("sparse.jl")

# TODO: use parameter to choose between Indef/PosDef system
AbstractLinearSolver(A::Matrix{Tv}) where{Tv<:Real} = DenseLinearSolver(A)
AbstractLinearSolver(
    A::SparseMatrixCSC{Tv, Int64}
) where{Tv<:BlasReal} = SparsePosDefLinearSolver(A)