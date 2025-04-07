# =============================================
#   1. Supported variable attributes
# =============================================

"""
    SUPPORTED_VARIABLE_ATTR

List of supported `MOI.VariableAttribute`.
* `MOI.VariablePrimal`
"""
const SUPPORTED_VARIABLE_ATTR = Union{
    MOI.VariableName,
    # MOI.VariablePrimalStart,
    MOI.VariablePrimal
}

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true


# =============================================
#   2. Add variables
# =============================================
function MOI.is_valid(m::Optimizer, x::MOI.VariableIndex)
    return haskey(m.var_indices, x)
end

function MOI.add_variable(m::Optimizer{T}) where{T}
    # TODO: dispatch a function call to m.inner instead of m.inner.pbdata
    m.var_counter += 1
    x = MOI.VariableIndex(m.var_counter)
    j = Tulip.add_variable!(m.inner.pbdata, Int[], T[], zero(T), T(-Inf), T(Inf))

    # Update tracking of variables
    m.var_indices[x] = j
    m.var2bndtype[x] = Set{Type{<:MOI.AbstractScalarSet}}()
    push!(m.var_indices_moi, x)

    return x
end

# =============================================
#   3. Delete variables
# =============================================
function MOI.delete(m::Optimizer, v::MOI.VariableIndex)
    MOI.throw_if_not_valid(m, v)

    # Update inner model
    j = m.var_indices[v]
    old_name = get_attribute(m.inner, VariableName(), j)
    delete_variable!(m.inner.pbdata, j)

    # Remove bound tracking
    delete!(m.var2bndtype, v)

    # Name update
    if old_name != ""
        s = m.name2var[old_name]
        delete!(s, v)
        length(s) == 0 && delete!(m.name2var, old_name)
    end

    # Update indices correspondence
    deleteat!(m.var_indices_moi, j)
    delete!(m.var_indices, v)
    for v_ in m.var_indices_moi[j:end]
        m.var_indices[v_] -= 1
    end
    return nothing
end

# =============================================
#   4. Get/set variable attributes
# =============================================

function MOI.get(m::Optimizer, ::Type{MOI.VariableIndex}, name::String)
    s = get(m.name2var, name, Set{MOI.VariableIndex}())
    if length(s) == 0
        return nothing
    elseif length(s) == 1
        return first(s)
    else
        error("Duplicate variable name detected: $(name)")
    end
end

function MOI.get(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    MOI.throw_if_not_valid(m, v)

    # Get name from inner model
    j = m.var_indices[v]
    return get_attribute(m.inner, VariableName(), j)
end

function MOI.set(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String)
    # Check that variable does exist
    MOI.throw_if_not_valid(m, v)

    # Update inner model
    j = m.var_indices[v]
    old_name = get_attribute(m.inner, VariableName(), j)
    if name == old_name
        return  # It's the same name!
    end
    set_attribute(m.inner, VariableName(), j, name)

    s = get!(m.name2var, name, Set{MOI.VariableIndex}())

    # Update names mapping
    push!(s, v)
    # Delete old name
    s_old = get(m.name2var, old_name, Set{MOI.ConstraintIndex}())
    if length(s_old) == 0
        # Variable didn't have name before
    elseif length(s_old) == 1
        # Delete this from mapping
        delete!(m.name2var, old_name)
    else
        delete!(s_old, v)
    end
    return nothing
end

function MOI.get(m::Optimizer{T},
    attr::MOI.VariablePrimal,
    x::MOI.VariableIndex
) where{T}
    MOI.throw_if_not_valid(m, x)
    MOI.check_result_index_bounds(m, attr)

    # Query inner solution
    j = m.var_indices[x]
    return m.inner.solution.x[j]
end
