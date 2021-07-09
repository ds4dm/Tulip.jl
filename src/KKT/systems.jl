"""
    DefaultKKTSystem

Default KKT system setting. Currently equivalent to [`K2`](@ref)
"""
struct DefaultKKTSystem <: AbstractKKTSystem end

@doc raw"""
    K2 <: AbstractKKTSystem

Augmented system
```math
    \begin{bmatrix}
        -(\Theta^{-1} + R_{p}) & A^{\top}\\
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
* ``\Theta^{-1} = X^{-l}Z^{l} + X^{-u} Z^{u}``
* ``R_{p}, R_{d}`` are current primal and dual regularizations
* ``\xi_{d}, \xi_{p}`` are given right-hand sides
"""
struct K2 <: AbstractKKTSystem end

@doc raw"""
    K1 <: AbstractKKTSystem

Normal equations system
```math
    \begin{array}{rl}
    \left(
        A (\Theta^{-1} + R_{p})^{-1} A^{\top} + R_{d}
    \right)
    \Delta y
    & =
    \xi_{p} + A (Θ^{-1} + R_{p})^{-1} \xi_{d}\\
    \Delta x &= (Θ^{-1} + R_{p})^{-1} (A^{\top} \Delta y - \xi_{d})
    \end{array}
```
where
* ``\Theta^{-1} = X^{-l}Z^{l} + X^{-u} Z^{u}``
* ``R_{p}, R_{d}`` are current primal and dual regularizations
* ``\xi_{d}, \xi_{p}`` are the augmented system's right-hand side
"""
struct K1 <: AbstractKKTSystem end
