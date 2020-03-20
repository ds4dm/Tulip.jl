# TODO: this is redundant with forcing constraint checks
#   => an empty row is automatically redundant or infeasible

struct EmptyRow{Tv} <: PresolveTransformation{Tv}
    i::Int  # row index
    y::Tv  # dual multiplier
end

function remove_empty_row!(lp::PresolveData{Tv}, i::Int) where{Tv}
    # Sanity checks
    (lp.rowflag[i] && lp.nzrow[i] == 0) || return nothing

    # Check bounds
    lb, ub = lp.lrow[i], lp.urow[i]

    if ub < zero(Tv)
        # Infeasible
        @info "Row $i is primal infeasible"
        lp.status = Trm_PrimalInfeasible
        lp.updated = true

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        resize!(lp.solution, lp.nrow, lp.ncol)
        lp.solution.x .= zero(Tv)
        lp.solution.y_lower .= zero(Tv)
        lp.solution.y_upper .= zero(Tv)
        lp.solution.s_lower .= zero(Tv)
        lp.solution.s_upper .= zero(Tv)

        lp.solution.primal_status = Sln_Unknown
        lp.solution.dual_status = Sln_InfeasibilityCertificate
        lp.solution.y_upper[i] = one(Tv)
        return
    elseif lb > zero(Tv)
        @info "Row $i is primal infeasible"
        lp.status = Trm_PrimalInfeasible
        lp.updated = true

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        resize!(lp.solution, lp.nrow, lp.ncol)
        lp.solution.x .= zero(Tv)
        lp.solution.y_lower .= zero(Tv)
        lp.solution.y_upper .= zero(Tv)
        lp.solution.s_lower .= zero(Tv)
        lp.solution.s_upper .= zero(Tv)
        
        lp.solution.primal_status = Sln_Unknown
        lp.solution.dual_status = Sln_InfeasibilityCertificate
        lp.solution.y_lower[i] = one(Tv)
        return
    else
        push!(lp.ops, EmptyRow(i, zero(Tv)))
    end

    # Book-keeping
    lp.updated = true
    lp.rowflag[i] = false
    lp.nrow -= 1
end

function postsolve!(sol::Solution{Tv}, sol_::Solution{Tv}, op::EmptyRow{Tv}) where{Tv}
    sol.y_lower[op.i] = pos_part(op.y)
    sol.y_upper[op.i] = neg_part(op.y)
    return nothing
end