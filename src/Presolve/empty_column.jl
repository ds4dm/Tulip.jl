struct EmptyColumn{T} <: PresolveTransformation{T}
    j::Int  # variable index
    x::T  # Variable value
    s::T  # Reduced cost
end

function remove_empty_column!(ps::PresolveData{T}, j::Int) where{T}
    # Sanity check
    (ps.colflag[j] && (ps.nzcol[j] == 0)) || return nothing

    # Remove column
    lb, ub = ps.lcol[j], ps.ucol[j]
    cj = ps.obj[j]
    @debug "Removing empty column $j" cj lb ub

    ϵ = ps.options.ToleranceDFeas

    if cj > ϵ
        if isfinite(lb)
            # Set variable to lower bound
            # Update objective constant
            ps.obj0 += lb * cj
            push!(ps.ops, EmptyColumn(j, lb, cj))
        else
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
    elseif cj < -ϵ
        if isfinite(ub)
            # Set variable to upper bound
            # Update objective constant
            ps.obj0 += ub * cj
            push!(ps.ops, EmptyColumn(j, ub, cj))
        else
            # Problem is dual infeasible
            @debug "Column $j is (upper) unbounded"
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

            # Unbounded ray: xj = 1
            ps.solution.primal_status = Sln_InfeasibilityCertificate
            ps.solution.dual_status = Sln_Unknown
            ps.solution.is_primal_ray = true
            ps.solution.is_dual_ray = false
            ps.solution.z_primal = ps.solution.z_dual = -T(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.x[j_] = one(T)

            return
        end
    else
        # Any feasible value works
        if isfinite(lb)
            push!(ps.ops, EmptyColumn(j, lb, zero(T)))
        elseif isfinite(ub)
            push!(ps.ops, EmptyColumn(j, ub, zero(T)))
        else
            # Free variable with zero coefficient and empty column
            push!(ps.ops, EmptyColumn(j, zero(T), zero(T)))
        end
    end

    # Book=keeping
    ps.colflag[j] = false
    ps.updated = true
    ps.ncol -= 1
    return nothing
end

function postsolve!(sol::Solution{T}, op::EmptyColumn{T}) where{T}
    sol.x[op.j] = op.x
    sol.s_lower[op.j] = pos_part(op.s)
    sol.s_upper[op.j] = neg_part(op.s)
    return nothing
end
