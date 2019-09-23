# Tulip

| **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:----------------:|:----------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ds4dm.github.io/Tulip.jl/dev/) | [![Build Status](https://travis-ci.org/ds4dm/Tulip.jl.svg?branch=master)](https://travis-ci.org/ds4dm/Tulip.jl) | [![codecov.io](https://codecov.io/github/ds4dm/Tulip.jl/coverage.svg?branch=master)](http://codecov.io/github/ds4dm/Tulip.jl?branch=master)


## Overview
Tulip is an open-source interior-point solver for linear optimization, written in pure Julia.
It implements the homogeneous primal-dual interior-point algorithm with multiple centrality corrections, and therefore handles unbounded and infeasible problems.
Tulip’s main feature is that its algorithmic framework is disentangled from linear algebra implementations.
This allows to seamlessly integrate specialized routines for structured problems.

## Installation

Just install like any Julia package

```julia
] add Tulip
```

## Usage

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
