"""
    ConstrId

"""
struct ConstrId
    uuid::Int64  # Unique identifier
end


"""
    LinConstrData{Tv<:Real}

"""
mutable struct LinConstrData{Tv<:Real}
    name::String  # Constraint name

    bt::BoundType

    lb::Tv
    ub::Tv

    function LinConstrData{Tv}(name::String, lb::Tv, ub::Tv) where{Tv<:Real}
        bt = bound_type(lb, ub)
        return new{Tv}(name, bt, lb, ub)
    end
end

function LinConstrData{Tv}(name::String, lb::Real, ub::Real) where{Tv<:Real}
    return LinConstrData{Tv}(name, Tv(lb), Tv(ub))
end

abstract type AbstractConstraint{Tv<:Real} end


"""
    LinearConstraint{Tv<:Real}

Place-holder for a linear constraint.

A linear constraint writes ``l \\leq  a^{T} x \\leq u``.
"""
struct LinearConstraint{Tv<:Real} <: AbstractConstraint{Tv}
    id::ConstrId      # Unique identifier
    dat::LinConstrData{Tv}  # Constraint data

    # Constructor
    # TODO: check
    function LinearConstraint(id::ConstrId, dat::LinConstrData{Tv}) where{Tv<:Real}
        return new{Tv}(id, dat)
    end
end

function LinearConstraint{Tv}(id::ConstrId, name::String, lb::Tv, ub::Tv) where{Tv<:Real}
    cdat = LinConstrData{Tv}(name, lb, ub)
    return LinearConstraint(id, cdat)
end

function LinearConstraint{Tv}(id::ConstrId, name::String, lb::Real, ub::Real) where{Tv<:Real}
    return LinearConstraint{Tv}(id, name, Tv(lb), Tv(ub))
end


"""
    get_uuid(c::LinearConstraint)

Return the identifier of variable `v`.
"""
get_uuid(c::LinearConstraint) = c.id


"""
    get_name(c::LinearConstraint)

Return the name of variable `v`.
"""
get_name(c::LinearConstraint) = c.dat.name


"""
    set_name(c::LinearConstraint, s::String)

Set the name of variable `v` to `s`.
"""
set_name!(c::LinearConstraint, s::String) = (c.dat.name = s)


"""
    get_bounds(c::LinearConstraint{Tv})
"""
function get_bounds(c::LinearConstraint{Tv}) where{Tv<:Real}
    return (c.dat.bt, c.dat.lb, c.dat.ub)
end


"""
    get_lower_bound(c::LinearConstraint)

"""
get_lower_bound(c::LinearConstraint) = c.dat.lb


"""
    get_upper_bound(c::LinearConstraint)

"""
get_upper_bound(c::LinearConstraint) = c.dat.ub


"""
    set_bounds!(c::LinearConstraint, bt, lb, ub)

"""
function set_bounds!(c::LinearConstraint{Tv}, lb::Tv, ub::Tv) where{Tv<:Real}
    # Check bounds
    bt = bound_type(lb, ub)

    c.dat.bt = bt
    c.dat.lb = lb
    c.dat.ub = ub
    return nothing
end

function set_bounds!(c::LinearConstraint{Tv}, lb::Real, ub::Real) where{Tv<:Real}
    set_bounds!(c, Tv(lb), Tv(ub))
end