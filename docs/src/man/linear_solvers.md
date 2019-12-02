```@meta
CurrentModule = Tulip.TLPLinearAlgebra
```

# Linear solvers

## Linear systems

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


## Overview

Here is a list of currently supported linear solvers:

| Linear solver type | type of ``A`` | System | Method | 
|:--------------------|:------:|:--:|:--|
| `DenseLinearSolver` | `Matrix` | Normal Eqn | Direct (Cholesky) |
| `SparseIndefLinearSolver` | `SparseMatricCSC` | Augm. Sys | Direct (LDLt) |
| `SparsePosDefLinearSolver` | `SparseMatricCSC` | Normal Eqn | Direct (Cholesky) |

## Linear solvers

```@docs
AbstractLinearSolver
```

Custom linear solvers should inherit from the `AbstractLinearSolver` class,
and extend the following two functions

```@docs
update_linear_solver!(::AbstractLinearSolver, ::Any, ::Any, ::Any)
```

```@docs
solve_augmented_system!(::Any, ::Any, ::AbstractLinearSolver, ::Any, ::Any)
```

### Dense

```@docs
DenseLinearSolver
```

```@docs
update_linear_solver!(::DenseLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:BlasReal}
```


### Sparse

```@docs
SparseIndefLinearSolver
```

```@docs
SparsePosDefLinearSolver
```