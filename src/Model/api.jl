"""
    add_constraint!(m::Model_)

Add one linear constraint to the model.
"""
function add_constraint!(m::Model_{Tv},
    name::String,
    lb::Tv, ub::Tv,
    vinds::Vector{VarId},
    vvals::Vector{Tv}
) where{Tv<:Real}

    # Sanity checks
    length(vinds) == length(vvals) || throw(DimensionMismatch(
        "vidx has length $(length(vinds)) but vval has length $(length(vvals))"
    ))

    # Check that all variables do exist
    for vidx in vinds
        haskey(m.pbdata_raw.vars, vidx) || error("Variable $(vidx.uuid) not in model.")
    end

    # Create new constraint object
    cidx = new_constraint_index!(m.pbdata_raw)
    constr = LinearConstraint{Tv}(cidx, name, lb, ub)

    # Add constraint to model
    add_constraint!(m.pbdata_raw, constr)

    # Set coefficients
    for (vidx, val) in zip(vinds, vvals)
        m.pbdata_raw.coeffs[vidx, cidx] = val
        push!(m.pbdata_raw.var2con[vidx], cidx)
        push!(m.pbdata_raw.con2var[cidx], vidx)
    end
    
    return cidx
end


"""
    add_variable!(m::Model_, name, obj, lb, ub, cinds, cvals)

Add one scalar variable to the model.
"""
function add_variable!(
    m::Model_{Tv},
    name::String,
    obj::Tv, lb::Tv, ub::Tv,
    cinds::Vector{ConstrId},
    cvals::Vector{Tv},
) where{Tv<:Real}

    # Sanity checks
    length(cinds) == length(cvals) || throw(DimensionMismatch(
        "vidx has length $(length(cinds)) but vval has length $(length(cvals))"
    ))

    # Check that all variables do exist
    for cidx in cinds
        haskey(m.pbdata_raw.constrs, cidx) || error("Constraint $(cidx.uuid) not in model.")
    end

    # Create new constraint object
    vidx = new_variable_index!(m.pbdata_raw)
    var = Variable{Tv}(vidx, name, obj, lb, ub)

    # Add constraint to model
    add_variable!(m.pbdata_raw, var)

    # Set coefficients
    for (cidx, val) in zip(cinds, cvals)
        m.pbdata_raw.coeffs[vidx, cidx] = val
        push!(m.pbdata_raw.var2con[vidx], cidx)
        push!(m.pbdata_raw.con2var[cidx], vidx)
    end
    
    return vidx
end

function add_variable!(
    m::Model_{Tv},
    name::String,
    obj::Tv, lb::Tv, ub::Tv
) where{Tv<:Real}

    # Create new Variable object
    vidx = new_variable_index!(m.pbdata_raw)
    var = Variable{Tv}(vidx, name, obj, lb, ub)

    # Add variable to model
    add_variable!(m.pbdata_raw, var)

    # No coefficient to add
    
    return vidx
end

function add_variable!(
    m::Model_{Tv}
) where{Tv<:Real}

    # Create new Variable object
    vidx = new_variable_index!(m.pbdata_raw)
    var = Variable{Tv}(vidx)

    # Add variable to model
    add_variable!(m.pbdata_raw, var)

    # No coefficient to add
    
    return vidx
end


"""
    add_constraints!(m::Model_{Tv}, ncons)

Add several linear constraints to the model.
"""
# TODO

"""
    set_constr_name!(m::Model_, i::Ti, name::String)

Set the name of a linear constraint.
"""
# function set_constr_name!(m::Model_, i::Ti, name::String) end


"""
    delete_constr!(m::Model_, i::Ti)

Delete a constraint from the model.
"""
# function delete_constr!(m::Model_, i::Ti) end


"""
    add_variable!(m::Model_)

Add one variable to the model.
"""
# function add_variable!(m::Model_) end


"""
    set_var_name!(m::Model_, j::Ti, name::String)

Set the name of a linear constraint.
"""
# function set_var_name!(m::Model_, j::Ti, name::String) end


"""
    delete_var!(m::Model_, j::Ti)

Delete a variable from the model.
"""
# function delete_var!(m::Model_, j::Ti) end


"""
    set_coefficient!(m::Model_, i::Ti, j::Ti, v::Tv)

Set the value of a matrix coefficient.
"""
# TODO: value of index types for Model_
# function set_coefficient!(m::Model_, i::Ti, j::Ti, v::Tv) end


"""
    set_constr_lb!(m::Model_, i::Ti, v::Tv)

Set the lower bound on a linear constraint.
"""
# function set_constr_lb!(m::Model_, i::Ti, v::Tv) end


"""
    set_constr_ub!(m::Model_, i::Ti, v::Tv)

Set the upper bound on a linear constraint.
"""
# function set_constr_ub!(m::Model_, i::Ti, v::Tv) end


"""
    set_var_lb!(m::Model_, i::Ti, v::Tv)

Set the lower bound on a scalar variable.
"""
# function set_var_lb!(m::Model_, j::Ti, v::Tv) end


"""
    set_var_ub!(m::Model_, i::Ti, v::Tv)

Set the upper bound on a scalar variable.
"""
# function set_var_ub!(m::Model_, j::Ti, v::Tv) end


"""
    set_objective_sense!(m::Model_, sense::)

Set the objective sense for the problem.
"""
# function set_objective_sense!(m::Model_, sense::)