"""
    VarId

"""
struct VarId
    uuid::Int64  # Unique identifier
end


"""
    VarData{Tv<:Real}

"""
mutable struct VarData{Tv<:Real}
    name::String

    obj::Tv  # Objective coeff

    bt::BoundType
    lb::Tv   # Lower bound
    ub::Tv   # Upper bound

    # Constructors
    function VarData{Tv}(name::String, obj::Real, bt::BoundType, lb::Real, ub::Real) where {Tv<:Real}
        return new{Tv}(name, obj, bt, lb, ub)
    end

end

function VarData(name::String, obj::Tv, bt::BoundType, lb::Tv, ub::Tv) where{Tv<:Real}
    VarData{Tv}(name, obj, bt, lb, ub)
end


"""
    Variable{Tv<:Real}

Place-holder for variables.
"""
struct Variable{Tv<:Real}
    id::VarId     # Unique identifier
    dat::VarData{Tv}  # Variable data

    # Constructor
    function Variable(id::VarId, dat::VarData{Tv}) where{Tv<:Real}
        return new{Tv}(id, dat)
    end
end

function Variable{Tv}(id::VarId, name::String, obj, bt, lb, ub) where{Tv<:Real}
    _check_bounds(bt, lb, ub) || error("Invalid bounds for $bt: [$lb, $ub]")
    vd = VarData{Tv}(name, obj, bt, lb, ub)
    return Variable(id, vd)
end


"""
    get_uuid(v::Variable)

Return the identifier of variable `v`.
"""
get_uuid(v::Variable) = v.id


"""
    get_name(v::Variable)

Return the name of variable `v`.
"""
get_name(v::Variable) = v.dat.name


"""
    set_name!(v::Variable, s::String)

Set the name of variable `v` to `s`.
"""
set_name!(v::Variable, s::String) = (v.dat.name = s)


"""
    get_obj_coeff(v::Variable)

Return the objective coefficient of variable `v`.
"""
get_obj_coeff(v::Variable) = v.dat.obj


"""
    set_obj_coeff(v::Variable, coeff)

Set objective coefficient of variable `v` to `c`.
"""
set_obj_coeff!(v::Variable, coeff) = (v.dat.obj = coeff)


"""
    get_bounds(v::Variable{Tv})
"""
function get_bounds(v::Variable{Tv}) where{Tv<:Real}
    return (v.dat.bt, v.dat.lb, v.dat.ub)
end

"""
    get_lower_bound(v::Variable)

"""
get_lower_bound(v::Variable) = v.dat.lb


"""
    get_upper_bound(v::Variable)

"""
get_upper_bound(v::Variable) = v.dat.ub


"""
    set_bounds!(c::LinearConstraint, bt, lb, ub)

"""
function set_bounds!(v::Variable{Tv}, bt::BoundType, lb::Tv, ub::Tv) where{Tv<:Real}
    # Check bounds
    _check_bounds(bt, lb, ub) || error("Invalid bounds for $bt: [$lb, $ub]")

    v.dat.bt = bt
    v.dat.lb = lb
    v.dat.ub = ub

    return nothing
end

set_bounds!(v::Variable{Tv}, bt, lb::Real, ub::Real) where {Tv<:Real}= set_bounds!(v, bt, Tv(lb), Tv(ub))