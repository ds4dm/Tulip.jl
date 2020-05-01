```@meta
CurrentModule = Tulip.KKT
```

# Solving linear systems

The interior-point algorithm in Tulip requires the solution, at each iteration, of the following symmetric quasi-definite _augmented system_
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
where
* ``\Delta x, \Delta y`` are primal and dual search directions,
* ``A`` is the problem's constraint matrix,
* ``\Theta``, ``R_p`` and ``R_d`` are positive diagonal matrices,
* ``\xi_p, \xi_d`` are right-hand side vectors.

The augmented system above can be reduced to the positive-definite _normal equations system_
```math
\begin{array}{rl}
\left(
    A (\Theta^{-1} + R_{p})^{-1} A^{T} + R_{d}
\right)
\Delta y
& =
\xi_{p} + A (Θ^{-1} + R_{p})^{-1} \xi_{d}\\
\Delta x &= (Θ^{-1} + R_{p})^{-1} (A^{T} \Delta y - \xi_{d})
\end{array}
```
When selected, this reduction is transparent to the interior-point algorithm.


To enable the use of fast external libraries and/or specialized routines, the resolution of linear systems is performed by a _linear solver_ object.
Linear solvers can be customized with the options below.

## Linear solver options

### System

```@docs
LinearSystem
```

```@docs
DefaultSystem
```

```@docs
AugmentedSystem
```

```@docs
NormalEquations
```


### Backend

```@docs
LSBackend
```

```@docs
DefaultBackend
```

```@docs
Lapack
```

```@docs
Cholmod
```

```@docs
LDLFact
```

### 


## Supported linear solvers

Here is a list of currently supported linear solvers:

| Linear solver type | `Tv` | System | `LSBackend` | Method |
|:-------------------|:----:|:------:|:-------:|:-------|
| [`KKTSolver_Dense`](@ref) | `Real` | `NormalEquations` | `Lapack` | Cholesky
| [`KKTSolver_CholmodQD`](@ref) | `Float64` | `AugmentedSystem` | `Cholmod` | LDLᵀ
| [`KKTSolver_CholmodPD`](@ref) | `Float64` | `NormalEquations` | `Cholmod` | Cholesky
| [`KKTSolver_LDLFact`](@ref) | `Real` | `AugmentedSystem` | `LDLFact` | LDLᵀ

### Default options
If no option is specified, then the linear solver is chosen as follows:
* ``A`` is dense: `KKTSolver_Dense`
* If ``A`` sparse and `Tv` is `Float64`: `KKTSolver_CholmodQD`
* All other cases: `KKTSolver_LDLFact`