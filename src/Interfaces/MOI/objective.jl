# =============================================
#   1. Supported objectives
# =============================================
MOI.supports(::Optimizer{Tv}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}) where{Tv} = true

# =============================================
#   2. Get/set objective function
# =============================================
function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}
) where{Tv}
    # Objective coeffs
    terms = MOI.ScalarAffineTerm{Tv}[]
    for (j, cj) in enumerate(m.inner.pbdata.obj)
        !iszero(cj) && push!(terms, MOI.ScalarAffineTerm(cj, m.var_indices_moi[j]))
    end

    # Constant term
    c0 = m.inner.pbdata.obj0

    return MOI.ScalarAffineFunction(terms, c0)
end

# TODO: use inner API
function MOI.set(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}},
    f::MOI.ScalarAffineFunction{Tv}
) where{Tv}

    # Sanity checks
    isfinite(f.constant) || error("Objective constant term must be finite")
    for t in f.terms
        MOI.throw_if_not_valid(m, t.variable_index)
    end

    # Update inner model
    m.inner.pbdata.obj .= zero(Tv) # Reset inner objective to zero
    for t in f.terms
        j = m.var_indices[t.variable_index]
        m.inner.pbdata.obj[j] += t.coefficient  # there may be dupplicates
    end
    set_attribute(m.inner, ObjectiveConstant(), f.constant)  # objective offset

    return nothing
end

# =============================================
#   3. Modify objective
# ============================================= 
function MOI.modify(
    m::Optimizer{Tv},
    c::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}},
    chg::MOI.ScalarCoefficientChange{Tv}
) where{Tv}
    # Sanity checks
    v = chg.variable
    MOI.throw_if_not_valid(m, v)

    # Update inner model
    j = m.var_indices[v]
    m.inner.pbdata.obj[j] = chg.new_coefficient  # TODO: use inner API
    return nothing
end

function MOI.modify(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}},
    chg::MOI.ScalarConstantChange{Tv}
) where{Tv}
    isfinite(chg.new_constant) || error("Objective constant term must be finite")
    m.inner.pbdata.obj0 = chg.new_constant
    return nothing
end