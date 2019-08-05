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

    lb::Tv
    ub::Tv

    # Q: add type (Equal, Less, Greater, Range)?
    LinConstrData{Tv}() where{Tv<:Real} = new{Tv}("", typemin(Tv), typemax(Tv))

    function LinConstrData{Tv}(name::String, lb, ub) where{Tv<:Real}
        return new{Tv}(name, lb, ub)
    end
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
    
    function LinearConstraint{Tv}(id::ConstrId) where{Tv<:Real}
        return LinearConstraint(id, LinConstrData{Tv}())
    end
end

function LinearConstraint{Tv}(id::ConstrId, name::String, lb, ub) where{Tv<:Real}
    cdat = LinConstrData{Tv}(name, lb, ub)
    return LinearConstraint(id, cdat)
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
    get_lower_bound(c::LinearConstraint)

"""
get_lower_bound(c::LinearConstraint) = c.dat.lb


"""
    set_lower_bound(c::LinearConstraint, l)

"""
set_lower_bound!(c::LinearConstraint, l) = (c.dat.lb = l)


"""
    get_upper_bound(c::LinearConstraint)

"""
get_upper_bound(c::LinearConstraint) = c.dat.ub


"""
    set_upper_bound(c::LinearConstraint, u)

"""
set_upper_bound!(c::LinearConstraint, u) = (c.dat.ub = u)