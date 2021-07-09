```@meta
CurrentModule = Tulip.KKT
```

# KKT systems

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

