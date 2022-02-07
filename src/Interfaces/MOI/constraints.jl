# =============================================
#   1. Supported constraints and attributes
# =============================================

"""
    SUPPORTED_CONSTR_ATTR

List of supported MOI `ConstraintAttribute`.
"""
const SUPPORTED_CONSTR_ATTR = Union{
    MOI.ConstraintName,
    MOI.ConstraintPrimal,
    MOI.ConstraintDual,
    # MOI.ConstraintPrimalStart,
    # MOI.ConstraintDualStart, # once dual warm-start is supported
    # MOI.ConstraintBasisStatus,  # once cross-over is supported
    MOI.ConstraintFunction,
    MOI.ConstraintSet
}

MOI.supports(::Optimizer, ::A, ::Type{<:MOI.ConstraintIndex}) where{A<:SUPPORTED_CONSTR_ATTR} = true

# MOI boilerplate
function MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex{<:MOI.VariableIndex}})
    throw(MOI.VariableIndexConstraintNameError())
end

# Variable bounds
function MOI.supports_constraint(
    ::Optimizer{T}, ::Type{MOI.VariableIndex}, ::Type{S}
) where {T, S<:SCALAR_SETS{T}}
    return true
end

# Linear constraints
function MOI.supports_constraint(
    ::Optimizer{T}, ::Type{MOI.ScalarAffineFunction{T}}, ::Type{S}
) where {T, S<:SCALAR_SETS{T}}
    return true
end


function MOI.is_valid(
    m::Optimizer{T},
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where{T, S <:SCALAR_SETS{T}}
    v = MOI.VariableIndex(c.value)
    MOI.is_valid(m, v) || return false
    res = S ∈ m.var2bndtype[v]
    return res
end

function MOI.is_valid(
    m::Optimizer{T},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    return haskey(m.con_indices, c)
end

# =============================================
#   2. Add constraints
# =============================================

# TODO: make it clear that only finite bounds can be given in input.
# To relax variable bounds, one should delete the associated bound constraint.
function MOI.add_constraint(
    m::Optimizer{T},
    v::MOI.VariableIndex,
    s::MOI.LessThan{T}
) where {T}

    # Check that variable exists
    MOI.throw_if_not_valid(m, v)
    # Check if upper bound already exists
    if MOI.LessThan{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{T}, MOI.LessThan{T}}(v))
    elseif MOI.EqualTo{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{T}, MOI.LessThan{T}}(v))
    elseif MOI.Interval{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{T}, MOI.LessThan{T}}(v))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableUpperBound(), j, s.upper)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.LessThan{T})

    return MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{T}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{T},
    v::MOI.VariableIndex,
    s::MOI.GreaterThan{T}
) where{T}

    # Check that variable exists
    MOI.throw_if_not_valid(m, v)
    # Check if lower bound already exists
    if MOI.GreaterThan{T} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{T}, MOI.GreaterThan{T}}(v))
    elseif MOI.EqualTo{T} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{T}, MOI.GreaterThan{T}}(v))
    elseif MOI.Interval{T} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{T}, MOI.GreaterThan{T}}(v))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.lower)

    # Update upper-bound
    push!(m.var2bndtype[v], MOI.GreaterThan{T})

    return MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{T}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{T},
    v::MOI.VariableIndex,
    s::MOI.EqualTo{T}
) where{T}

    # Check that variable exists
    MOI.throw_if_not_valid(m, v)
    # Check if a bound already exists
    if MOI.LessThan{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{T}, MOI.EqualTo{T}}(v))
    elseif MOI.GreaterThan{T} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{T}, MOI.EqualTo{T}}(v))
    elseif MOI.EqualTo{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{T}, MOI.EqualTo{T}}(v))
    elseif MOI.Interval{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{T}, MOI.EqualTo{T}}(v))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.value)
    set_attribute(m.inner, VariableUpperBound(), j, s.value)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.EqualTo{T})

    return MOI.ConstraintIndex{MOI.VariableIndex, MOI.EqualTo{T}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{T},
    v::MOI.VariableIndex,
    s::MOI.Interval{T}
) where{T}

    # Check that variable exists
    MOI.throw_if_not_valid(m, v)
    # Check if a bound already exists
    if MOI.LessThan{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{T}, MOI.Interval{T}}(v))
    elseif MOI.GreaterThan{T} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{T}, MOI.Interval{T}}(v))
    elseif MOI.EqualTo{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{T}, MOI.Interval{T}}(v))
    elseif MOI.Interval{T} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{T}, MOI.Interval{T}}(v))
    end

    # Update variable bounds
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.lower)
    set_attribute(m.inner, VariableUpperBound(), j, s.upper)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.Interval{T})

    return MOI.ConstraintIndex{MOI.VariableIndex, MOI.Interval{T}}(v.value)
end

# General linear constraints
function MOI.add_constraint(
    m::Optimizer{T},
    f::MOI.ScalarAffineFunction{T},
    s::SCALAR_SETS{T}
) where{T}
    # Check that constant term is zero
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{T, typeof(f), typeof(s)}(f.constant))
    end

    # Convert to canonical form
    fc = MOI.Utilities.canonical(f)

    # Extract row
    nz = length(fc.terms)
    rind = Vector{Int}(undef, nz)
    rval = Vector{T}(undef, nz)
    lb, ub = _bounds(s)
    for (k, t) in enumerate(fc.terms)
        rind[k] = m.var_indices[t.variable]
        rval[k] = t.coefficient
    end

    # Update inner model
    i = add_constraint!(m.inner.pbdata, rind, rval, lb, ub)

    # Create MOI index
    m.con_counter += 1
    cidx = MOI.ConstraintIndex{typeof(f), typeof(s)}(m.con_counter)

    # Update constraint tracking
    m.con_indices[cidx] = i
    push!(m.con_indices_moi, cidx)

    return cidx
end

# =============================================
#   3. Delete constraints
# =============================================
function MOI.delete(
    m::Optimizer{T},
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where{T, S<:SCALAR_SETS{T}}

    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Update inner model
    j = m.var_indices[v]
    if S == MOI.LessThan{T}
        # Remove upper-bound
        set_attribute(m.inner, VariableUpperBound(), j, T(Inf))
    elseif S == MOI.GreaterThan{T}
        # Remove lower bound
        set_attribute(m.inner, VariableLowerBound(), j, T(-Inf))
    else
        # Set variable to free
        set_attribute(m.inner, VariableLowerBound(), j, T(-Inf))
        set_attribute(m.inner, VariableUpperBound(), j, T(Inf))
    end

    # Update name tracking
    old_name = get(m.bnd2name, c, "")
    if old_name != "" && haskey(m.name2con, old_name)
        s = m.name2con[old_name]
        delete!(s, c)
        length(s) == 0 && delete!(m.name2con, old_name)
    end
    delete!(m.bnd2name, c)

    # Delete tracking of bounds
    delete!(m.var2bndtype[v], S)

    return nothing
end

function MOI.delete(
    m::Optimizer{T},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)

    # Update inner model
    i = m.con_indices[c]
    old_name = get_attribute(m.inner, ConstraintName(), i)
    delete_constraint!(m.inner.pbdata, i)

    # Update index tracking
    for c_ in m.con_indices_moi[i+1:end]
        m.con_indices[c_] -= 1
    end
    deleteat!(m.con_indices_moi, i)
    delete!(m.con_indices, c)

    # Update name tracking
    if old_name != "" && haskey(m.name2con, old_name)
        s = m.name2con[old_name]
        delete!(s, c)
        length(s) == 0 && delete!(m.name2con, old_name)
    end

    return nothing
end

# =============================================
#   4. Modify constraints
# =============================================
function MOI.modify(
    m::Optimizer{T},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S},
    chg::MOI.ScalarCoefficientChange{T}
) where{T, S<:SCALAR_SETS{T}}
    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))
    MOI.is_valid(m, chg.variable) || throw(MOI.InvalidIndex(chg.variable))

    # Update inner problem
    i = m.con_indices[c]
    j = m.var_indices[chg.variable]
    v = chg.new_coefficient

    set_coefficient!(m.inner.pbdata, i, j, v)
    return nothing
end

# =============================================
#   5. Get/set constraint attributes
# =============================================

#
#   ListOfConstraintIndices
#
function MOI.get(
    m::Optimizer{T},
    ::MOI.ListOfConstraintIndices{MOI.VariableIndex, S}
) where{T, S<:SCALAR_SETS{T}}
    indices = MOI.ConstraintIndex{MOI.VariableIndex, S}[]

    for (var, bounds_set) in m.var2bndtype
        S ∈ bounds_set && push!(indices, MOI.ConstraintIndex{MOI.VariableIndex, S}(var.value))
    end
    return sort!(indices, by = v -> v.value)
end

function MOI.get(
    m::Optimizer{T},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    indices = [
        cidx
        for cidx in keys(m.con_indices) if isa(cidx,
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
        )
    ]
    return sort!(indices, by = v -> v.value)
end

#
#   NumberOfConstraints
#
function MOI.get(
    m::Optimizer{T},
    ::MOI.NumberOfConstraints{MOI.VariableIndex, S}
) where{T, S<:SCALAR_SETS{T}}
    ncon = 0
    for (v, bound_sets) in m.var2bndtype
        ncon += S ∈ bound_sets
    end
    return ncon
end

function MOI.get(
    m::Optimizer{T},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    ncon = 0

    for cidx in keys(m.con_indices)
        ncon += isa(cidx, MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S})
    end

    return ncon
end

#
#   ConstraintName
#

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where {T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)

    return get(m.bnd2name, c, "")
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where {T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)

    # Get name from inner model
    i = m.con_indices[c]
    return get_attribute(m.inner, ConstraintName(), i)
end

function MOI.set(::Optimizer, ::MOI.ConstraintName, ::MOI.ConstraintIndex{<:MOI.VariableIndex}, ::String)
    throw(MOI.VariableIndexConstraintNameError())
end

function MOI.set(
    m::Optimizer{T}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S},
    name::String
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)

    s = get!(m.name2con, name, Set{MOI.ConstraintIndex}())

    # Update inner model
    i = m.con_indices[c]
    old_name = get_attribute(m.inner, ConstraintName(), i)
    set_attribute(m.inner, ConstraintName(), i, name)

    # Update constraint name tracking
    push!(s, c)
    # Delete old name
    s_old = get(m.name2con, old_name, Set{MOI.ConstraintIndex}())
    if length(s_old) == 0
        # Constraint previously didn't have name --> ignore
    elseif length(s_old) == 1
        delete!(m.name2con, old_name)
    else
        delete!(s_old, c)
    end
    return nothing
end

function MOI.get(m::Optimizer, CIType::Type{<:MOI.ConstraintIndex}, name::String)
    s = get(m.name2con, name, Set{MOI.ConstraintIndex}())
    if length(s) == 0
        return nothing
    elseif length(s) == 1
        c = first(s)
        return isa(c, CIType) ? c : nothing
    else
        error("Duplicate constraint name detected: $(name)")
    end
    return nothing
end

#
#   ConstraintFunction
#
function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where {T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    return MOI.VariableIndex(c.value)
end

function MOI.set(
    ::Optimizer{T}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S},
    ::MOI.VariableIndex,
) where {T, S<:SCALAR_SETS{T}}
    return throw(MOI.SettingVariableIndexNotAllowed())
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where {T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get row from inner model
    i = m.con_indices[c]
    row = m.inner.pbdata.arows[i]
    nz = length(row.nzind)

    # Map inner indices to MOI indices
    terms = Vector{MOI.ScalarAffineTerm{T}}(undef, nz)
    for (k, (j, v)) in enumerate(zip(row.nzind, row.nzval))
        terms[k] = MOI.ScalarAffineTerm{T}(v, m.var_indices_moi[j])
    end

    return MOI.ScalarAffineFunction(terms, zero(T))
end

# TODO
function MOI.set(
    m::Optimizer{T}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S},
    f::MOI.ScalarAffineFunction{T}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)
    iszero(f.constant) || throw(MOI.ScalarFunctionConstantNotZero{T, typeof(f), S}(f.constant))

    fc = MOI.Utilities.canonical(f)

    # Update inner model
    # TODO: use inner query
    i = m.con_indices[c]
    # Set old row to zero
    f_old = MOI.get(m, MOI.ConstraintFunction(), c)
    for term in f_old.terms
        j = m.var_indices[term.variable]
        set_coefficient!(m.inner.pbdata, i, j, zero(T))
    end
    # Set new row coefficients
    for term in fc.terms
        j = m.var_indices[term.variable]
        set_coefficient!(m.inner.pbdata, i, j, term.coefficient)
    end

    # Done

    return nothing
end

#
#   ConstraintSet
#
function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{T}}
) where{T}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    ub = m.inner.pbdata.uvar[j]

    return MOI.LessThan(ub)
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{T}}
) where{T}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    lb = m.inner.pbdata.lvar[j]

    return MOI.GreaterThan(lb)
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex, MOI.EqualTo{T}}
) where{T}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    ub = m.inner.pbdata.uvar[j]

    return MOI.EqualTo(ub)
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex, MOI.Interval{T}}
) where{T}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    lb = m.inner.pbdata.lvar[j]
    ub = m.inner.pbdata.uvar[j]

    return MOI.Interval(lb, ub)
end

function MOI.get(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get inner bounds
    i = m.con_indices[c]
    lb = m.inner.pbdata.lcon[i]
    ub = m.inner.pbdata.ucon[i]

    if S == MOI.LessThan{T}
        return MOI.LessThan(ub)
    elseif S == MOI.GreaterThan{T}
        return MOI.GreaterThan(lb)
    elseif S == MOI.EqualTo{T}
        return MOI.EqualTo(lb)
    elseif S == MOI.Interval{T}
        return MOI.Interval(lb, ub)
    end
end

function MOI.set(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S},
    s::S
) where{T, S<:SCALAR_SETS{T}}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Update inner bounds
    # Bound key does not need to be updated
    j = m.var_indices[v]
    if S == MOI.LessThan{T}
        set_attribute(m.inner, VariableUpperBound(), j, s.upper)
    elseif S == MOI.GreaterThan{T}
        set_attribute(m.inner, VariableLowerBound(), j, s.lower)
    elseif S == MOI.EqualTo{T}
        set_attribute(m.inner, VariableLowerBound(), j, s.value)
        set_attribute(m.inner, VariableUpperBound(), j, s.value)
    elseif S == MOI.Interval{T}
        set_attribute(m.inner, VariableLowerBound(), j, s.lower)
        set_attribute(m.inner, VariableUpperBound(), j, s.upper)
    else
        error("Unknown type for ConstraintSet: $S.")
    end

    return nothing
end

function MOI.set(
    m::Optimizer{T}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S},
    s::S
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)

    # Update inner bounds
    i = m.con_indices[c]
    if S == MOI.LessThan{T}
        set_attribute(m.inner, ConstraintUpperBound(), i, s.upper)
    elseif S == MOI.GreaterThan{T}
        set_attribute(m.inner, ConstraintLowerBound(), i, s.lower)
    elseif S == MOI.EqualTo{T}
        set_attribute(m.inner, ConstraintLowerBound(), i, s.value)
        set_attribute(m.inner, ConstraintUpperBound(), i, s.value)
    elseif S == MOI.Interval{T}
        set_attribute(m.inner, ConstraintLowerBound(), i, s.lower)
        set_attribute(m.inner, ConstraintUpperBound(), i, s.upper)
    else
        error("Unknown type for ConstraintSet: $S.")
    end

    return nothing
end

#
#   ConstraintPrimal
#
function MOI.get(
    m::Optimizer{T}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Query row primal
    j = m.var_indices[MOI.VariableIndex(c.value)]
    return m.inner.solution.x[j]
end

function MOI.get(
    m::Optimizer{T}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Query from inner model
    i = m.con_indices[c]
    return m.inner.solution.Ax[i]
end

#
#   ConstraintDual
#
function MOI.get(
    m::Optimizer{T}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.VariableIndex, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Get variable index
    j = m.var_indices[MOI.VariableIndex(c.value)]

    # Extract reduced cost
    if S == MOI.LessThan{T}
        return -m.inner.solution.s_upper[j]
    elseif S == MOI.GreaterThan{T}
        return m.inner.solution.s_lower[j]
    else
        return m.inner.solution.s_lower[j] - m.inner.solution.s_upper[j]
    end
end

function MOI.get(
    m::Optimizer{T}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, S}
) where{T, S<:SCALAR_SETS{T}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Get dual from inner model
    i = m.con_indices[c]
    if isa(S, MOI.LessThan)
        return -m.inner.solution.y_upper[i]
    elseif isa(S, MOI.GreaterThan)
        return m.inner.solution.y_lower[i]
    else
        return m.inner.solution.y_lower[i] - m.inner.solution.y_upper[i]
    end
end
