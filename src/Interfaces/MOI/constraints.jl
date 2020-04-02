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

# Variable bounds
function MOI.supports_constraint(
    ::Optimizer{Tv}, ::Type{MOI.SingleVariable}, ::Type{S}
) where {Tv, S<:SCALAR_SETS{Tv}}
    return true
end

# Linear constraints
function MOI.supports_constraint(
    ::Optimizer{Tv}, ::Type{MOI.ScalarAffineFunction{Tv}}, ::Type{S}
) where {Tv, S<:SCALAR_SETS{Tv}}
    return true
end


function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S <:SCALAR_SETS{Tv}}
    v = MOI.VariableIndex(c.value)
    MOI.is_valid(m, v) || return false
    res = S ∈ m.var2bndtype[v]
    return res
end

function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    return haskey(m.con_indices, c)
end

# =============================================
#   2. Add constraints
# =============================================

# TODO: make it clear that only finite bounds can be given in input.
# To relax variable bounds, one should delete the associated bound constraint.
function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.LessThan{Tv}
) where{Tv}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    # Check if upper bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.LessThan{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.LessThan{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.LessThan{Tv}}(v.value))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableUpperBound(), j, s.upper)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.LessThan{Tv})
    
    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.GreaterThan{Tv}
) where{Tv}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    # Check if lower bound already exists
    if MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.GreaterThan{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Tv}, MOI.GreaterThan{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Tv}, MOI.GreaterThan{Tv}}(v.value))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.lower)

    # Update upper-bound
    push!(m.var2bndtype[v], MOI.GreaterThan{Tv})

    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.EqualTo{Tv}
) where{Tv}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    # Check if a bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.EqualTo{Tv}}(v.value))
    end

    # Update inner model
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.value)
    set_attribute(m.inner, VariableUpperBound(), j, s.value)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.EqualTo{Tv})

    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}(v.value)
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.Interval{Tv}
) where{Tv}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    # Check if a bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.Interval{Tv}}(v.value))
    end

    # Update variable bounds
    j = m.var_indices[v]  # inner index
    set_attribute(m.inner, VariableLowerBound(), j, s.lower)
    set_attribute(m.inner, VariableUpperBound(), j, s.upper)

    # Update bound tracking
    push!(m.var2bndtype[v], MOI.Interval{Tv})

    return MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}(v.value)
end

# General linear constraints
function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.ScalarAffineFunction{Tv},
    s::S
) where{Tv, S<:SCALAR_SETS{Tv}}
    # Check that constant term is zero
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Tv, typeof(f), S}(f.constant))
    end

    # Convert to canonical form
    fc = MOI.Utilities.canonical(f)

    # Extract row
    nz = length(fc.terms)
    rind = Vector{Int}(undef, nz)
    rval = Vector{Tv}(undef, nz)
    lb, ub = _bounds(s)
    for (k, t) in enumerate(fc.terms)
        rind[k] = m.var_indices[t.variable_index]
        rval[k] = t.coefficient
    end

    # Update inner model
    i = add_constraint!(m.inner.pbdata, rind, rval, lb, ub)

    # Create MOI index
    m.con_counter += 1
    cidx = MOI.ConstraintIndex{typeof(f), S}(m.con_counter)

    # Update constraint tracking
    m.con_indices[cidx] = i
    push!(m.con_indices_moi, cidx)

    return cidx
end

# =============================================
#   3. Delete constraints
# =============================================
function MOI.delete(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}

    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Update inner model
    j = m.var_indices[v]
    if S == MOI.LessThan{Tv}
        # Remove upper-bound
        set_attribute(m.inner, VariableUpperBound(), j, Tv(Inf))
    elseif S == MOI.GreaterThan{Tv}
        # Remove lower bound
        set_attribute(m.inner, VariableLowerBound(), j, Tv(-Inf))
    else
        # Set variable to free
        set_attribute(m.inner, VariableLowerBound(), j, Tv(-Inf))
        set_attribute(m.inner, VariableUpperBound(), j, Tv(Inf))
    end

    # Update name tracking
    old_name = get(m.bnd2name, c, "")
    old_name != "" && delete!(m.name2con, old_name)  # delete old name
    delete!(m.bnd2name, c)

    # Delete tracking of bounds
    delete!(m.var2bndtype[v], S)

    return nothing
end

function MOI.delete(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
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
    old_name != "" && delete!(m.name2con, old_name)

    return nothing
end

# =============================================
#   4. Modify constraints
# =============================================
function MOI.modify(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    chg::MOI.ScalarCoefficientChange{Tv}
) where{Tv, S<:SCALAR_SETS{Tv}}
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
    m::Optimizer{Tv},
    ::MOI.ListOfConstraintIndices{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[]

    for (var, bounds_set) in m.var2bndtype
        S ∈ bounds_set && push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(var.value))
    end
    return indices
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    return [
        cidx
        for cidx in keys(m.con_indices) if isa(cidx,
            MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
        )
    ]
end

#
#   NumberOfConstraints
#
function MOI.get(
    m::Optimizer{Tv},
    ::MOI.NumberOfConstraints{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    ncon = 0
    for (v, bound_sets) in m.var2bndtype
        ncon += S ∈ bound_sets
    end
    return ncon
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    ncon = 0

    for cidx in keys(m.con_indices)
        ncon += isa(cidx, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S})
    end
    
    return ncon
end

#
#   ConstraintName
#

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    
    return get(m.bnd2name, c, "")
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    # Get name from inner model
    i = m.con_indices[c]
    return get_attribute(m.inner, ConstraintName(), i)
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    name::String
) where{Tv, S<:SCALAR_SETS{Tv}}
    # Sanity checks
    MOI.throw_if_not_valid(m, c)
    c_ = get(m.name2con, name, nothing)
    c_ === nothing || c_ == c || error("Dupplicate constraint name $name")

    # Remove old name
    old_name = get(m.bnd2name, c, "")
    delete!(m.name2con, old_name)

    # Update new name
    if name == ""
        delete!(m.bnd2name, c)
    else
        m.bnd2name[c] = name
        m.name2con[name] = c
    end
    return nothing
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    name::String
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    # Check for dupplicate name
    c_ = get(m.name2con, name, nothing)
    c_ === nothing || c_ == c || error("Dupplicate constraint name $name")

    # Update inner model
    i = m.con_indices[c]
    old_name = get_attribute(m.inner, ConstraintName(), i)
    set_attribute(m.inner, ConstraintName(), i, name)

    # Update constraint name tracking
    delete!(m.name2con, old_name)
    if name != ""
        m.name2con[name] = c
    end
    return nothing
end

function MOI.get(m::Optimizer, CIType::Type{<:MOI.ConstraintIndex}, name::String)
    c = get(m.name2con, name, nothing)
    return isa(c, CIType) ? c : nothing
end

#
#   ConstraintFunction
#
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    ::MOI.SingleVariable
) where{Tv, S<:SCALAR_SETS{Tv}}
    return throw(MOI.SettingSingleVariableFunctionNotAllowed())
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get row from inner model
    i = m.con_indices[c]
    row = m.inner.pbdata.arows[i]
    nz = length(row.nzind)

    # Map inner indices to MOI indices
    terms = Vector{MOI.ScalarAffineTerm{Tv}}(undef, nz)
    for (k, (j, v)) in enumerate(zip(row.nzind, row.nzval))
        terms[k] = MOI.ScalarAffineTerm{Tv}(v, m.var_indices_moi[j])
    end

    return MOI.ScalarAffineFunction(terms, zero(Tv))
end

# TODO
function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    f::MOI.ScalarAffineFunction{Tv}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    iszero(f.constant) || throw(MOI.ScalarFunctionConstantNotZero{Tv, typeof(f), S}(f.constant))

    fc = MOI.Utilities.canonical(f)

    # Update inner model
    # TODO: use inner query
    i = m.con_indices[c]
    # Set old row to zero
    f_old = MOI.get(m, MOI.ConstraintFunction(), c)
    for term in f_old.terms
        j = m.var_indices[term.variable_index]
        set_coefficient!(m.inner.pbdata, i, j, zero(Tv))
    end
    # Set new row coefficients
    for term in fc.terms
        j = m.var_indices[term.variable_index]
        set_coefficient!(m.inner.pbdata, i, j, term.coefficient)
    end

    # Done

    return nothing
end

#
#   ConstraintSet
#
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}
) where{Tv}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    ub = m.inner.pbdata.uvar[j]

    return MOI.LessThan(ub)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}
) where{Tv}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    lb = m.inner.pbdata.lvar[j]

    return MOI.GreaterThan(lb)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}
) where{Tv}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Get inner bounds
    j  = m.var_indices[v]
    ub = m.inner.pbdata.uvar[j]

    return MOI.EqualTo(ub)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}
) where{Tv}
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
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get inner bounds
    i = m.con_indices[c]
    lb = m.inner.pbdata.lcon[i]
    ub = m.inner.pbdata.ucon[i]

    if S == MOI.LessThan{Tv}
        return MOI.LessThan(ub)
    elseif S == MOI.GreaterThan{Tv}
        return MOI.GreaterThan(lb)
    elseif S == MOI.EqualTo{Tv}
        return MOI.EqualTo(lb)
    elseif S == MOI.Interval{Tv}
        return MOI.Interval(lb, ub)
    end
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    s::S
) where{Tv, S<:SCALAR_SETS{Tv}}
    # Sanity check
    MOI.throw_if_not_valid(m, c)
    v = MOI.VariableIndex(c.value)

    # Update inner bounds
    # Bound key does not need to be updated
    j = m.var_indices[v]
    if S == MOI.LessThan{Tv}
        set_attribute(m.inner, VariableUpperBound(), j, s.upper)
    elseif S == MOI.GreaterThan{Tv}
        set_attribute(m.inner, VariableLowerBound(), j, s.lower)
    elseif S == MOI.EqualTo{Tv}
        set_attribute(m.inner, VariableLowerBound(), j, s.value)
        set_attribute(m.inner, VariableUpperBound(), j, s.value)
    elseif S == MOI.Interval{Tv}
        set_attribute(m.inner, VariableLowerBound(), j, s.lower)
        set_attribute(m.inner, VariableUpperBound(), j, s.upper)
    else
        error("Unknown type for ConstraintSet: $S.")
    end

    return nothing
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    s::S
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    
    # Update inner bounds
    i = m.con_indices[c]
    if S == MOI.LessThan{Tv}
        set_attribute(m.inner, ConstraintUpperBound(), i, s.upper)
    elseif S == MOI.GreaterThan{Tv}
        set_attribute(m.inner, ConstraintLowerBound(), i, s.lower)
    elseif S == MOI.EqualTo{Tv}
        set_attribute(m.inner, ConstraintLowerBound(), i, s.value)
        set_attribute(m.inner, ConstraintUpperBound(), i, s.value)
    elseif S == MOI.Interval{Tv}
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
    m::Optimizer{Tv}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Query row primal
    j = m.var_indices[MOI.VariableIndex(c.value)]
    return m.inner.solution.col_primal[j]
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Query from inner model
    i = m.con_indices[c]
    return m.inner.solution.row_primal[i]
end

#
#   ConstraintDual
#
function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Get variable index
    j = m.var_indices[MOI.VariableIndex(c.value)]
    rc = m.inner.solution.col_dual[j]

    # Extract reduced cost
    if S == MOI.LessThan{Tv}
        return min(rc, 0)
    elseif S == MOI.GreaterThan{Tv}
        return max(rc, 0)
    else
        return rc
    end
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    # Get dual from inner model
    i = m.con_indices[c]
    return m.inner.solution.row_dual[i]
end