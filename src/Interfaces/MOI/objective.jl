# =============================================
#   1. Supported objectives
# =============================================
function MOI.supports(
    ::Optimizer{T},
    ::MOI.ObjectiveFunction{F}
) where{T, F<:Union{MOI.VariableIndex, MOI.ScalarAffineFunction{T}}}
    return true
end

# =============================================
#   2. Get/set objective function
# =============================================
function MOI.get(
    m::Optimizer{T},
    ::MOI.ObjectiveFunction{F}
) where{T,F}
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
    return convert(F, obj)
end

function MOI.get(
    m::Optimizer{T},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}
) where{T}
    # Objective coeffs
    terms = MOI.ScalarAffineTerm{T}[]
    for (j, cj) in enumerate(m.inner.pbdata.obj)
        !iszero(cj) && push!(terms, MOI.ScalarAffineTerm(cj, m.var_indices_moi[j]))
    end

    # Constant term
    c0 = m.inner.pbdata.obj0

    return MOI.ScalarAffineFunction(terms, c0)
end

# TODO: use inner API
function MOI.set(
    m::Optimizer{T},
    ::MOI.ObjectiveFunction{F},
    f::F
) where{T, F <: MOI.VariableIndex}

    MOI.set(
        m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(),
        convert(MOI.ScalarAffineFunction{T}, f)
    )
    m._obj_type = _SINGLE_VARIABLE
    return nothing
end

function MOI.set(
    m::Optimizer{T},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    f::MOI.ScalarAffineFunction{T}
) where{T}

    # Sanity checks
    isfinite(f.constant) || error("Objective constant term must be finite")
    for t in f.terms
        MOI.throw_if_not_valid(m, t.variable)
    end

    # Update inner model
    m.inner.pbdata.obj .= zero(T) # Reset inner objective to zero
    for t in f.terms
        j = m.var_indices[t.variable]
        m.inner.pbdata.obj[j] += t.coefficient  # there may be dupplicates
    end
    set_attribute(m.inner, ObjectiveConstant(), f.constant)  # objective offset

    m._obj_type = _SCALAR_AFFINE

    return nothing
end

# =============================================
#   3. Modify objective
# =============================================
function MOI.modify(
    m::Optimizer{T},
    c::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    chg::MOI.ScalarCoefficientChange{T}
) where{T}
    # Sanity checks
    v = chg.variable
    MOI.throw_if_not_valid(m, v)

    # Update inner model
    j = m.var_indices[v]
    m.inner.pbdata.obj[j] = chg.new_coefficient  # TODO: use inner API
    m._obj_type = _SCALAR_AFFINE
    return nothing
end

function MOI.modify(
    m::Optimizer{T},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    chg::MOI.ScalarConstantChange{T}
) where{T}
    isfinite(chg.new_constant) || error("Objective constant term must be finite")
    m.inner.pbdata.obj0 = chg.new_constant
    m._obj_type = _SCALAR_AFFINE
    return nothing
end
