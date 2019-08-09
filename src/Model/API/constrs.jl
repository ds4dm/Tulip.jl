# ==================================
#       ADD LINEAR CONSTRAINTS
# ==================================

"""
    add_constraint!(m::Model_)

Add one linear constraint to the model.
"""
function add_constraint!(
    m::Model_{Tv},
    name::String,
    bt::BoundType, lb::Tv, ub::Tv,
    colids::Vector{VarId},
    colvals::Vector{Tv}
) where{Tv<:Real}

    # =================
    #   Sanity checks
    # =================
    _check_bounds(bt, lb, ub) || error("Invalid bounds for $bt: [$lb, $ub].")

    length(colids) == length(colvals) || throw(DimensionMismatch(
        "colidx has length $(length(colids)) but vval has length $(length(colvals))"
    ))

    # Check that all variables do exist
    for colidx in colids
        haskey(m.pbdata_raw.vars, colidx) || error("Variable $(colidx.uuid) not in model.")
    end


    # ============================
    #   Create Constraint object
    # ============================
    cidx = new_constraint_index!(m.pbdata_raw)
    constr = LinearConstraint{Tv}(cidx, name, bt, lb, ub)

    # Add constraint to model
    add_constraint!(m.pbdata_raw, constr)


    # ====================
    #   Set coefficients
    # ====================
    # TODO: put this in a function
    # In particular, explicit zeros might be erasing an existing non-zero
    for (colidx, val) in zip(colids, colvals)
        set_coeff!(m.pbdata_raw, colidx, cidx, val)
    end
    
    return cidx
end

add_constraint!(
    m::Model_{Tv},
    name::String,
    bt::BoundType, lb::Real, ub::Real,
    colids::Vector{VarId},
    colvals::Vector{T}
) where{Tv<:Real, T<:Real} = add_constraint!(
    m, name,
    bt, Tv(lb), Tv(ub),
    colids, Tv.(colvals)
)

add_constraint!(
    m::Model_{Tv},
    name::String,
    obj::Real,
    bt::BoundType, lb::Real, ub::Real
) where{Tv<:Real} = add_constraint!(m, name, obj, bt, lb, ub, VarrId[], Tv[])