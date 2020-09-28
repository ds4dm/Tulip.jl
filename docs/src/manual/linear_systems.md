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
If available and when selected, this reduction is transparent to the interior-point algorithm.

To enable the use of fast external libraries and/or specialized routines, the resolution of linear systems is performed by an [`AbstractKKTSolver`] object.


## Supported linear solvers

Here is a list of currently supported linear solvers:

| Linear solver type | System | Backend | Method |
|:-------------------|:-------|:--------|:-------|
| [`DenseSPD`](@ref) | Normal equations | Dense / LAPACK | Cholesky
| [`CholmodSQD`](@ref) | Augmented system | CHOLMOD | LDLᵀ
| [`CholmodSPD`](@ref) | Normal equations | CHOLMOD | Cholesky
| [`LDLFactSQD`](@ref) | Augmented system | [LDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) | LDLᵀ
| [`KrylovSPD`](@ref) | Normal equations | [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) | Krylov
| [`KrylovSID`](@ref) | Augmented system[^1] | [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) | Krylov
| [`KrylovSQD`](@ref) | Augmented system[^1] | [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) | Krylov

[^1]: [`KrylovSID`](@ref)s view the augmented system as a symmetric indefinite system,
    while [`KrylovSQD`](@ref)s exploit its 2x2 structure and quasi-definite property.
    See the reference documentation for more details.