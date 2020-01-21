```@meta
CurrentModule = Tulip.TLPLinearAlgebra
```

# Solving linear systems

!!! warning

    The overall structure of Linear Solvers is expected to undergo a major update in the near future.
    In particular, a number of type parametrization shall eventually be replaced by a trait-based implementation, 
    which should prove more flexible.

    The high-level interface should not change too much, and the underlying linear algebra
        techniques shall remain the same.

## Overview

The first system we consider is the symmetric quasi-definite _augmented system_
```math
\begin{bmatrix}
    -(\Theta^{-1} + R_{p}) & A^{T}\\
    A & R_{d}
\end{bmatrix}
\begin{bmatrix}
    \Delta x\\
    \Delta y
\end{bmatrix}
=
\begin{bmatrix}
    \xi_d\\
    \xi_p
\end{bmatrix}
```

One can pivot out the upper-left diagonal block to obtain the positive-definite
    _normal equations system_
```math
\left(
    A (\Theta^{-1} + R_{p})^{-1} A^{T} + R_{d}
\right)
\Delta y
=
\xi_{p} + A (Θ^{-1} + R_{p})^{-1} \xi_{d}
```


## Linear solvers

Here is a list of currently supported linear solvers:

| Linear solver type | type of ``A`` | System | Method | Backend |
|:--------------------|:------:|:--:|:--|:--|
| [`DenseLinearSolver`](@ref) | `Matrix` | Normal Eqn | Direct (Cholesky) | LAPACK |
| [`SparseIndefLinearSolver`](@ref) | `SparseMatricCSC` | Augm. Sys | Direct (LDLt) | SuiteSparse |
| [`SparsePosDefLinearSolver`](@ref) | `SparseMatricCSC` | Normal Eqn | Direct (Cholesky) | SuiteSparse |
| [`LDLFLinearSolver`](@ref) | `SparseMatricCSC` | Augm. Sys | Direct (LDLt) | [LDLFactorizations](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) |

### AbstractLinearSolver

This is the base type from which all implementations should derive.

```@docs
AbstractLinearSolver
```

Custom linear solvers should inherit from the `AbstractLinearSolver` class,
and extend the following two functions:

```@docs
update_linear_solver!
```

```@docs
solve_augmented_system!
```

### DenseLinearSolver

```@docs
DenseLinearSolver
```

```@docs
update_linear_solver!(::DenseLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
update_linear_solver!(::DenseLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:BlasReal}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::DenseLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::DenseLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:BlasReal}
```

### SparseIndefLinearSolver

```@docs
SparseIndefLinearSolver
```

```@docs
update_linear_solver!(::SparseIndefLinearSolver,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::SparseIndefLinearSolver, ::Vector{Float64}, ::Vector{Float64})
```



### SparsePosDefLinearSolver

```@docs
SparsePosDefLinearSolver
```

```@docs
update_linear_solver!(::SparsePosDefLinearSolver,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::SparsePosDefLinearSolver, ::Vector{Float64}, ::Vector{Float64})
```

### LDLFLinearSolver

```@docs
LDLFLinearSolver
```

```@docs
update_linear_solver!(::LDLFLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::LDLFLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
```