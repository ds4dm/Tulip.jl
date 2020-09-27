# TODO: this is redundant with forcing constraint checks
#   => an empty row is automatically redundant or infeasible

struct EmptyRow{T} <: PresolveTransformation{T}
    i::Int  # row index
    y::T  # dual multiplier
end

function remove_empty_row!(ps::PresolveData{T}, i::Int) where{T}
    # Sanity checks
    (ps.rowflag[i] && ps.nzrow[i] == 0) || return nothing

    # Check bounds
    lb, ub = ps.lrow[i], ps.urow[i]

    if ub < zero(T)
        # Infeasible
        @debug "Row $i is primal infeasible"
        ps.status = Trm_PrimalInfeasible
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(T)
        ps.solution.y_lower .= zero(T)
        ps.solution.y_upper .= zero(T)
        ps.solution.s_lower .= zero(T)
        ps.solution.s_upper .= zero(T)

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = Sln_Unknown
        ps.solution.dual_status = Sln_InfeasibilityCertificate
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = T(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_upper[i] = one(T)
        return
    elseif lb > zero(T)
        @debug "Row $i is primal infeasible"
        ps.status = Trm_PrimalInfeasible
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(T)
        ps.solution.y_lower .= zero(T)
        ps.solution.y_upper .= zero(T)
        ps.solution.s_lower .= zero(T)
        ps.solution.s_upper .= zero(T)
        
        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = Sln_Unknown
        ps.solution.dual_status = Sln_InfeasibilityCertificate
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = T(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_lower[i_] = one(T)
        return
    else
        push!(ps.ops, EmptyRow(i, zero(T)))
    end

    # Book-keeping
    ps.updated = true
    ps.rowflag[i] = false
    ps.nrow -= 1
end

function postsolve!(sol::Solution{T}, op::EmptyRow{T}) where{T}
    sol.y_lower[op.i] = pos_part(op.y)
    sol.y_upper[op.i] = neg_part(op.y)
    return nothing
end