```@meta
CurrentModule = Tulip.KKT
```

# Overview

The `KKT` module provides a modular, customizable interface for developing and selecting various approaches to solve the KKT systems.

## KKT backends

```@docs
AbstractKKTBackend
```

```@docs
DefaultKKTBackend
```

## KKT systems

All formulations below refer to a linear program in primal-dual standard form
```math
    \begin{array}{rl}
    (P) \ \ \ 
    \displaystyle \min_{x}
    & c^{\top} x\\
    s.t.
    & A x = b\\
    & l \leq x \leq u
    \end{array}
    \quad \quad \quad
    \begin{array}{rl}
    (D) \ \ \ 
    \displaystyle \max_{y, z}
    & b^{\top} y + l^{\top}z^{l} - u^{\top}z^{u}\\
    s.t.
    & A^{\top}y + z^{l} - z^{u} = c\\
    & z^{l}, z^{u} \geq 0
    \end{array}
```

```@docs
AbstractKKTSystem
```

```@docs
DefaultKKTSystem
```

```@docs
K2
```

```@docs
K1
```

## KKT solvers

```@docs
AbstractKKTSolver
```

Custom linear solvers should (preferably) inherit from the `AbstractKKTSolver` class,
and extend the following functions:

```@docs
setup
```

```@docs
update!
```

```@docs
solve!
```