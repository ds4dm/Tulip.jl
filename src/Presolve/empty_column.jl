struct EmptyColumn{Tv} <: PresolveTransformation{Tv}
    j::Int  # variable index
    x::Tv  # Variable value
    s::Tv  # Reduced cost
end

function remove_empty_column!(lp::PresolveData{Tv}, j::Int) where{Tv}
    # Sanity check
    (lp.colflag[j] && (lp.nzcol[j] == 0)) || return nothing

    # Remove column
    lb, ub = lp.lcol[j], lp.ucol[j]
    cj = lp.obj[j]
    @debug "Removing empty column $j" cj lb ub

    if cj > zero(Tv)
        if isfinite(lb)
            # Set variable to lower bound
            # Update objective constant
            lp.obj0 += lb * cj
            push!(lp.ops, EmptyColumn(j, lb, cj))
        else
            # Problem is dual infeasible
            @debug "Column $j is (lower) unbounded"
            lp.status = Trm_DualInfeasible
            lp.updated = true

            # Resize problem
            compute_index_mapping!(lp)
            resize!(lp.solution, lp.nrow, lp.ncol)
            lp.solution.x .= zero(Tv)
            lp.solution.y_lower .= zero(Tv)
            lp.solution.y_upper .= zero(Tv)
            lp.solution.s_lower .= zero(Tv)
            lp.solution.s_upper .= zero(Tv)

            # Unbounded ray: xj = -1
            lp.solution.primal_status = Sln_InfeasibilityCertificate
            lp.solution.dual_status = Sln_Unknown
            lp.solution.z_primal = lp.solution.z_dual = -Tv(Inf)
            j_ = lp.new_var_idx[j]
            lp.solution.x[j_] = -one(Tv)

            return nothing
        end
    elseif cj < zero(Tv)
        if isfinite(ub)
            # Set variable to upper bound
            # Update objective constant
            lp.obj0 += ub * cj
            push!(lp.ops, EmptyColumn(j, ub, cj))
        else
            # Problem is dual infeasible
            @debug "Column $j is (upper) unbounded"
            lp.status = Trm_DualInfeasible
            lp.updated = true

            # Resize problem
            compute_index_mapping!(lp)
            resize!(lp.solution, lp.nrow, lp.ncol)
            lp.solution.x .= zero(Tv)
            lp.solution.y_lower .= zero(Tv)
            lp.solution.y_upper .= zero(Tv)
            lp.solution.s_lower .= zero(Tv)
            lp.solution.s_upper .= zero(Tv)

            # Unbounded ray: xj = 1
            lp.solution.primal_status = Sln_InfeasibilityCertificate
            lp.solution.dual_status = Sln_Unknown
            lp.solution.z_primal = lp.solution.z_dual = -Tv(Inf)
            j_ = lp.new_var_idx[j]
            lp.solution.x[j_] = one(Tv)
            
            return
        end
    else
        # Any feasible value works
        if isfinite(lb)
            push!(lp.ops, EmptyColumn(j, lb, zero(Tv)))
        elseif isfinite(ub)
            push!(lp.ops, EmptyColumn(j, ub, zero(Tv)))
        else
            # Free variable with zero coefficient and empty column
            push!(lp.ops, EmptyColumn(j, zero(Tv), zero(Tv)))
        end

    end

    # Book=keeping
    lp.colflag[j] = false
    lp.updated = true
    lp.ncol -= 1
    return nothing
end

function postsolve!(sol::Solution{Tv}, sol_::Solution{Tv}, op::EmptyColumn{Tv}) where{Tv}
    sol.x[op.j] = op.x
    sol.s_lower[op.j] = pos_part(op.s)
    sol.s_upper[op.j] = neg_part(op.s)
    return nothing
end