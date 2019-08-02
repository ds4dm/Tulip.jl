"""
    Variable{Tv<:Real, Ti<:Integer}

Place-holder for variables.
"""
mutable struct Variable{Tv<:Real, Ti<:Integer}
    uuid::Ti      # Unique identifier
    name::String  # Variable name

    obj::Tv  # Objective coeff

    lb::Tv  # Lower bound
    ub::Tv  # Upper bound

    # Constructor
    # TODO: check and add useful constructors
    function Variable{Tv, Ti}(uuid, name, obj, lb, ub) where{Tv<:Real, Ti<:Integer}
        return new{Tv, Ti}(uuid, name, obj, lb, ub)
    end
end


"""
    Variable{Tv, Ti}(uuid::Integer)

"""
function Variable{Tv, Ti}(uuid::Integer) where{Tv<:Real, Ti<:Integer}
    return Variable{Tv, Ti}(uuid, "", zero(Tv), typemin(Tv), typemax(Tv))
end


"""
    get_uuid(v::Variable)

Return the identifier of variable `v`.
"""
get_uuid(v::Variable) = v.uuid


"""
    get_name(v::Variable)

Return the name of variable `v`.
"""
get_name(v::Variable) = v.name


"""
    set_name(v::Variable, s::String)

Set the name of variable `v` to `s`.
"""
set_name!(v::Variable, s::String) = (v.name = s)


"""
    get_obj_coeff(v::Variable)

Return the objective coefficient of variable `v`.
"""
get_obj_coeff(v::Variable) = v.obj


"""
    set_obj_coeff(v::Variable, coeff)

Set objective coefficient of variable `v` to `c`.
"""
set_obj_coeff!(v::Variable, coeff) = (v.obj = coeff)


"""
    get_lower_bound(v::Variable)

"""
get_lower_bound(v::Variable) = v.lb


"""
    set_lower_bound(v::Variable, l)

"""
set_lower_bound!(v::Variable, l) = (v.lb = l)


"""
    get_upper_bound(v::Variable)

"""
get_upper_bound(v::Variable) = v.ub


"""
    set_upper_bound(v::Variable, u)

"""
set_upper_bound!(v::Variable, u) = (v.ub = u)


"""
    get_rows(v::Variable)

"""
get_rows(v::Variable) = copy(v.rows)


"""
    add_row(v::Variable, r)

Add one row to variable `v`.
"""
add_row!(v::Variable, r) = push!(v.rows, r)


"""
    delete_row!(v::Variable, r)

Delete row `r` from variable `v`.
"""
delete_row!(v::Variable, r) = error("delete_row! not implemented yet.")
