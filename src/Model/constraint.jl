abstract type AbstractConstraint{Tv<:Real, Ti<:Integer} end


"""
    LinearConstraint{Tv<:Real, Ti<:Integer}

Place-holder for a linear constraint.

A linear constraint writes ``l \\leq  a^{T} x \\leq u``.
"""
mutable struct LinearConstraint{Tv<:Real, Ti<:Integer} <: AbstractConstraint{Tv, Ti}
    uuid::Ti      # Unique identifier
    name::String  # Constraint name

    lb::Tv  # Lower bound
    ub::Tv  # Upper bound
    # Q: add type (Equal, Less, Greater, Range)?

    # Indices of variables that appear in this constraint
    # Actual coefficients are stored in a separate container
    # Q: Do these indices have to be sorted?
    # Q: Do we even need those indices here?
    cols::OrderedSet{Ti}

    # Constructor
    # TODO: check
    function LinearConstraint{Tv, Ti}(uuid, name, lb, ub, cols) where{Tv<:Real, Ti<:Integer}
        return new{Tv, Ti}(uuid, name, lb, ub, cols)
    end
end