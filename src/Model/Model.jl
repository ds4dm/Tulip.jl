"""
    Point{T, Tv}

Data structure for primal-dual point.

Numerical precision is `T`, while `Tv` is the vector type.
"""
mutable struct Point{T<:Real, Tv<:AbstractVector{T}}
        
    # Primal-dual coordinates
    x::Tv  # primal variables
    w::Tv  # primal upper-bound slack

    y::Tv  # dual variables
    s::Tv  # reduced cost of `x`
    z::Tv  # reduced cost of `w`

    t::T  # Homogeneous primal variable
    k::T  # Homogeneous dual variable

end

"""
    Model{Ta, Tv, Ti}

Place-holder for problem data and related info.
"""
mutable struct Model{Ta<:AbstractMatrix{<:Real}, Tv<:AbstractVector{<:Real}, Ti<:AbstractVector{<:Integer}}
    A::Ta  # Constraints matrix
    b::Tv  # Right-hand side
    c::Tv  # Objective
    uind::Ti  # Indices of upper-bounded variables
    uval::Tv  # Finite upper bounds on variables
end

