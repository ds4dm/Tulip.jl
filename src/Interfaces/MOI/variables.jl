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

# TODO: dispatch to inner model
function MOI.add_variables(m::Optimizer, N::Int)
    N >= 0 || error("Cannot add negative number of variables")

    N == 0 && return MOI.VariableIndex[]

    vars = Vector{MOI.VariableIndex}(undef, N)
    for j in 1:N
        x = MOI.add_variable(m)
        vars[j] = x
    end

    return vars
end


# =============================================
#   3. Delete variables
# =============================================
function MOI.delete(m::Optimizer, v::MOI.VariableIndex)
    MOI.throw_if_not_valid(m, v)

    # Update inner model
    j = m.var_indices[v]
    delete_variable!(m.inner.pbdata, j)
    
    # Remove bound tracking
    delete!(m.var2bndtype, v)
    
    # TODO: name update

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

MOI.get(m::Optimizer, ::Type{MOI.VariableIndex}, name::String) = get(m.name2var, name, nothing)

function MOI.get(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    MOI.throw_if_not_valid(m, v)

    # Get name from inner model
    j = m.var_indices[v]
    return get_attribute(m.inner, VariableName(), j)
end

function MOI.set(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String)
    # Check that variable does exist
    MOI.throw_if_not_valid(m, v)

    # Check that name is unique
    v_ = get(m.name2var, name, nothing)
    v_ === nothing || v_ == v || error("Duplicate variable name $name")

    # Update inner model
    j = m.var_indices[v]
    set_attribute(m.inner, VariableName(), j, name)

    # Update names mapping
    m.name2var[name] = v
    
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