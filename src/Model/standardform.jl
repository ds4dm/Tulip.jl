using SparseArrays

"""
    StandardForm{Tv, Ti, Ta}

Problem data converted in standard form.
```math
\\begin{align}
    \\min_{x} \\ \\ \\ & c^{T} x \\\\
    s.t. \\ \\ \\
    & A x = b \\\\
    & x_{i} \\leq u \\\\
    & x \\geq 0
\\end{align}
```
"""
mutable struct StandardForm{Tv<:Real}
    A::AbstractMatrix{Tv}  # Constraint matrix
    b::Vector{Tv}          # Right-hand side
    c::Vector{Tv}          # Objective
    uind::Vector{Int}      # Indices of upper-bounded variables
    uval::Vector{Tv}       # Finite upper bounds on variables

    # TODO: add optimization sense
    # TODO: add row and column scalings
    # TODO: add starting points
    # TODO: add Cholesky factor (?)
end