struct DominatedColumn{Tv} <: PresolveTransformation{Tv}
    j::Int
    x::Tv  # Primal value
    cj::Tv  # Objective
    col::Col{Tv}  # Column
end

function remove_dominated_column!(lp::PresolveData{Tv}, j::Int; tol::Tv=100*sqrt(eps(Tv))) where{Tv}
    lp.colflag[j] || return nothing

    # Compute implied bounds on reduced cost: `ls ≤ s ≤ us`
    ls = us = zero(Tv)
    col = lp.pb0.acols[j]
    for (i, aij) in zip(col.nzind, col.nzval)
        (lp.rowflag[i] && !iszero(aij)) || continue

        ls += aij * ( (aij >= zero(Tv)) ? lp.ly[i] : lp.uy[i] )
        us += aij * ( (aij >= zero(Tv)) ? lp.uy[i] : lp.ly[i] )
    end

    # Check if column is dominated
    cj = lp.obj[j]
    if cj - us > tol
        # Reduced cost is always positive => fix to lower bound (or problem is unbounded)
        lb = lp.lcol[j]
        @debug "Fixing dominated column $j to its lower bound $lb"

        if !isfinite(lb)
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
            lp.solution.is_primal_ray = true
            lp.solution.is_dual_ray = false
            lp.solution.z_primal = lp.solution.z_dual = -Tv(Inf)
            j_ = lp.new_var_idx[j]
            lp.solution.x[j_] = -one(Tv)

            return nothing
        end

        # Update objective
        lp.obj0 += cj * lb

        # Extract column and update rows
        col_ = Col{Tv}(Int[], Tv[])
        for (i, aij) in zip(col.nzind, col.nzval)
            lp.rowflag[i] || continue

            push!(col_.nzind, i)
            push!(col_.nzval, aij)

            # Update bounds and non-zeros
            lp.lrow[i] -= aij * lb
            lp.urow[i] -= aij * lb
            lp.nzrow[i] -= 1

            lp.nzrow[i] == 1 && push!(lp.row_singletons, i)
        end

        # Remove variable from problem
        push!(lp.ops, DominatedColumn(j, lb, cj, col_))
        lp.colflag[j] = false
        lp.ncol -= 1
        lp.updated = true

    elseif cj - ls < -tol
        # Reduced cost is always negative => fix to upper bound (or problem is unbounded)
        ub = lp.ucol[j]
        
        if !isfinite(ub)
            # Problem is unbounded
            @debug "Column $j is (upper) unbounded"

            lp.status = Trm_DualInfeasible
            lp.updated = true

            # Resize solution
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
            lp.solution.is_primal_ray = true
            lp.solution.is_dual_ray = false
            lp.solution.z_primal = lp.solution.z_dual = -Tv(Inf)
            j_ = lp.new_var_idx[j]
            lp.solution.x[j_] = one(Tv)

            return nothing
        end

        @debug "Fixing dominated column $j to its upper bound $ub"

        # Update objective
        lp.obj0 += cj * ub

        # Extract column and update rows
        col_ = Col{Tv}(Int[], Tv[])
        for (i, aij) in zip(col.nzind, col.nzval)
            lp.rowflag[i] || continue

            push!(col_.nzind, i)
            push!(col_.nzval, aij)

            # Update bounds and non-zeros
            lp.lrow[i] -= aij * ub
            lp.urow[i] -= aij * ub
            lp.nzrow[i] -= 1

            lp.nzrow[i] == 1 && push!(lp.row_singletons, i)
        end

        # Remove variable from problem
        push!(lp.ops, DominatedColumn(j, ub, cj, col_))
        lp.colflag[j] = false
        lp.ncol -= 1
        lp.updated = true
    end

    return nothing
end

function postsolve!(sol::Solution{Tv}, op::DominatedColumn{Tv}) where{Tv}
    # Primal value
    sol.x[op.j] = op.x

    # Reduced cost
    s = sol.is_dual_ray ? zero(Tv) : op.cj
    for (i, aij) in zip(op.col.nzind, op.col.nzval)
        s -= aij * (sol.y_lower[i] - sol.y_upper[i])
    end

    sol.s_lower[op.j] = pos_part(s)
    sol.s_upper[op.j] = neg_part(s)

    return nothing
end