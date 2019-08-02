"""
    StandardFormData{Tv, Ti, Ta}

Problem data converted in standard form.
"""
mutable struct StandardFormData{Tv<:Real, Ti<:Integer}
    A::AbstractMatrix{Tv}  # Constraint matrix
    b::Vector{Tv}          # Right-hand side
    c::Vector{Tv}          # Objective
    uind::Vector{Ti}       # Indices of upper-bounded variables
    uval::Vector{Tv}       # Finite upper bounds on variables

    # TODO: add optimization sense
    # TODO: add row and column scalings
    # TODO: add starting points
    # TODO: add Cholesky factor (?)
end

# TODO:
# Add function to convert from ProblemData to StandardFormData