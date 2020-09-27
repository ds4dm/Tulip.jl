# Problem formulation

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
    & l \leq x \leq u\\
    \end{array}
```
where ``x, c \in \mathbb{R}^{n}``, ``A \in \mathbb{R}^{m \times n}``, ``b \in \mathbb{R}^{m}``,
and  ``l, u \in (\mathbb{R} \cup \{-\infty, +\infty \})^{n}``, i.e., some bounds may be infinite.

The original problem is automatically reformulated into standard form before the optimization is performed.
This transformation is transparent to the user.