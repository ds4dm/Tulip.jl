# Formulations

## Model input

Tulip takes as input LP problems of the form
```math
    \begin{array}{rrcll}
    (P) \ \ \ 
    \displaystyle \min_{x} && c^{T}x & + \ c_{0}\\
    s.t.
    & l^{b}_{i} \leq & a_{i}^{T} x & \leq u^{b}_{i} & \forall i = 1, ..., m\\
    & l^{x}_{j} \leq & x_{j} & \leq u^{x}_{j} & \forall j = 1, ..., n\\
    \end{array}
```
where ``l^{b,x}, u^{b, x} \in \mathbb{R} \cup \{ - \infty, + \infty \}``, i.e., some of the bounds may be infinite.

This original formulation is then converted to standard form.

## Standard form

Internally, Tulip solves LPs of the form
```math
    \begin{array}{rl}
    (P) \ \ \ 
    \displaystyle \min_{x}
    & c^{T} x + \ c_{0}\\
    s.t.
    & A x = b\\
    & x \leq u\\
    & x \geq 0
    \end{array}
```
where ``x, c, u \in \mathbb{R}^{n}``, ``A \in \mathbb{R}^{m \times n}`` and ``b \in \mathbb{R}^{m}``.
Some ``u_{j}`` may may take infinite value, i.e., the corresponding variable ``x_{j}`` has no upper bound.

The original problem is automatically reformulated into standard form before the optimization is performed.
This transformation is transparent to the user.