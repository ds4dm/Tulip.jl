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
    lb::Tv   # Lower bound
    ub::Tv   # Upper bound

    # Constructors
    VarData{Tv}() where {Tv<:Real} = new{Tv}("", zero(Tv), typemin(Tv), typemax(Tv))
    
    function VarData{Tv}(name::String, obj, lb, ub) where {Tv<:Real}
        return new{Tv}(name, obj, lb, ub)
    end

end

VarData(name::String, obj::Tv, lb::Tv, ub::Tv) where{Tv<:Real} = VarData{Tv}(name, obj, lb, ub)


"""
    Variable{Tv<:Real}

Place-holder for variables.
"""
struct Variable{Tv<:Real}
    id::VarId     # Unique identifier
    dat::VarData{Tv}  # Variable data

    # Constructor
    function Variable(id::VarId, dat::VarData{Tv}) where{Tv<:Real, Ti<:Integer}
        return new{Tv}(id, dat)
    end
    Variable{Tv}(id::VarId) where{Tv<:Real, Ti<:Integer} = Variable(id, VarData{Tv}())
end

function Variable{Tv}(id::VarId, name::String, obj, lb, ub) where{Tv<:Real, Ti<:Integer}
    vd = VarData{Tv}(name, obj, lb, ub)
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
    set_name(v::Variable, s::String)

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
    get_lower_bound(v::Variable)

"""
get_lower_bound(v::Variable) = v.dat.lb


"""
    set_lower_bound(v::Variable, l)

"""
set_lower_bound!(v::Variable, l) = (v.dat.lb = l)


"""
    get_upper_bound(v::Variable)

"""
get_upper_bound(v::Variable) = v.dat.ub


"""
    set_upper_bound(v::Variable, u)

"""
set_upper_bound!(v::Variable, u) = (v.dat.ub = u)