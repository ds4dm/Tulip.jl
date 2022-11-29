struct DominatedColumn{T} <: PresolveTransformation{T}
    j::Int
    x::T  # Primal value
    cj::T  # Objective
    col::Col{T}  # Column
end

function remove_dominated_column!(ps::PresolveData{T}, j::Int; tol::T=100*sqrt(eps(T))) where{T}
    ps.colflag[j] || return nothing

    # Compute implied bounds on reduced cost: `ls ≤ s ≤ us`
    ls = us = zero(T)
    col = ps.pb0.acols[j]
    for (i, aij) in zip(col.nzind, col.nzval)
        (ps.rowflag[i] && !iszero(aij)) || continue

        ls += aij * ( !signbit(aij) ? ps.ly[i] : ps.uy[i] )
        us += aij * ( !signbit(aij) ? ps.uy[i] : ps.ly[i] )
    end

    # Check if column is dominated
    cj = ps.obj[j]
    if cj - us > tol
        # Reduced cost is always positive => fix to lower bound (or problem is unbounded)
        lb = ps.lcol[j]
        @debug "Fixing dominated column $j to its lower bound $lb"

        if !isfinite(lb)
            # Problem is dual infeasible
            @debug "Column $j is (lower) unbounded"
            ps.status = Trm_DualInfeasible
            ps.updated = true

            # Resize problem
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(T)
            ps.solution.y_lower .= zero(T)
            ps.solution.y_upper .= zero(T)
            ps.solution.s_lower .= zero(T)
            ps.solution.s_upper .= zero(T)

            # Unbounded ray: xj = -1
            ps.solution.primal_status = Sln_InfeasibilityCertificate
            ps.solution.dual_status = Sln_Unknown
            ps.solution.is_primal_ray = true
            ps.solution.is_dual_ray = false
            ps.solution.z_primal = ps.solution.z_dual = -T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.x[j_] = -one(T)

            return nothing
        end

        # Update objective
        ps.obj0 += cj * lb

        # Extract column and update rows
        col_ = Col{T}(Int[], T[])
        for (i, aij) in zip(col.nzind, col.nzval)
            ps.rowflag[i] || continue

            push!(col_.nzind, i)
            push!(col_.nzval, aij)

            # Update bounds and non-zeros
            ps.lrow[i] -= aij * lb
            ps.urow[i] -= aij * lb
            ps.nzrow[i] -= 1

            ps.nzrow[i] == 1 && push!(ps.row_singletons, i)
        end

        # Remove variable from problem
        push!(ps.ops, DominatedColumn(j, lb, cj, col_))
        ps.colflag[j] = false
        ps.ncol -= 1
        ps.updated = true

    elseif cj - ls < -tol
        # Reduced cost is always negative => fix to upper bound (or problem is unbounded)
        ub = ps.ucol[j]
        
        if !isfinite(ub)
            # Problem is unbounded
            @debug "Column $j is (upper) unbounded"

            ps.status = Trm_DualInfeasible
            ps.updated = true

            # Resize solution
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(T)
            ps.solution.y_lower .= zero(T)
            ps.solution.y_upper .= zero(T)
            ps.solution.s_lower .= zero(T)
            ps.solution.s_upper .= zero(T)

            # Unbounded ray: xj = -1
            ps.solution.primal_status = Sln_InfeasibilityCertificate
            ps.solution.dual_status = Sln_Unknown
            ps.solution.is_primal_ray = true
            ps.solution.is_dual_ray = false
            ps.solution.z_primal = ps.solution.z_dual = -T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.x[j_] = one(T)

            return nothing
        end

        @debug "Fixing dominated column $j to its upper bound $ub"

        # Update objective
        ps.obj0 += cj * ub

        # Extract column and update rows
        col_ = Col{T}(Int[], T[])
        for (i, aij) in zip(col.nzind, col.nzval)
            ps.rowflag[i] || continue

            push!(col_.nzind, i)
            push!(col_.nzval, aij)

            # Update bounds and non-zeros
            ps.lrow[i] -= aij * ub
            ps.urow[i] -= aij * ub
            ps.nzrow[i] -= 1

            ps.nzrow[i] == 1 && push!(ps.row_singletons, i)
        end

        # Remove variable from problem
        push!(ps.ops, DominatedColumn(j, ub, cj, col_))
        ps.colflag[j] = false
        ps.ncol -= 1
        ps.updated = true
    end

    return nothing
end

function postsolve!(sol::Solution{T}, op::DominatedColumn{T}) where{T}
    # Primal value
    sol.x[op.j] = op.x

    # Reduced cost
    s = sol.is_dual_ray ? zero(T) : op.cj
    for (i, aij) in zip(op.col.nzind, op.col.nzval)
        s -= aij * (sol.y_lower[i] - sol.y_upper[i])
    end

    sol.s_lower[op.j] = pos_part(s)
    sol.s_upper[op.j] = neg_part(s)

    return nothing
end
