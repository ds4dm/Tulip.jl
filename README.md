

# Tulip

Tulip is an Interior-Point solver  written entirely in Julia.
It uses a standard Primal-Dual Predictor-Corrector algorithm to solve Linear Programs (LPs).

[![Build Status](https://travis-ci.org/ds4dm/Tulip.jl.svg?branch=master)](https://travis-ci.org/ds4dm/Tulip.jl)

[![codecov.io](http://codecov.io/github/ds4dm/Tulip.jl/coverage.svg?branch=master)](http://codecov.io/github/ds4dm/Tulip.jl?branch=master)

# Overview

Tulip solves LPs of the form:

$$
    \begin{array}{rrl}
    (P) \ \ \ 
    \displaystyle \min_{x} & c^{T}x\\
    s.t.
    & Ax &=b\\
    & x & \leq u\\
    & x & \geq 0\\
    \end{array}
$$
where $x, c, u \in \mathbb{R}^{n}$, $A \in \mathbb{R}^{m \times n}$ and $b \in \mathbb{R}^{m}$.
Some $u_{i}$ may may take infinite value, i.e., the corresponding variable $x_{i}$ has no upper bound.

It uses a primal-dual predictor-corrector interior point.
