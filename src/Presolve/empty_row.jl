# TODO: this is redundant with forcing constraint checks
#   => an empty row is automatically redundant or infeasible

struct EmptyRow{Tv} <: PresolveTransformation{Tv}
    i::Int  # row index
    y::Tv  # dual multiplier
end

function remove_empty_row!(ps::PresolveData{Tv}, i::Int) where{Tv}
    # Sanity checks
    (ps.rowflag[i] && ps.nzrow[i] == 0) || return nothing

    # Check bounds
    lb, ub = ps.lrow[i], ps.urow[i]

    if ub < zero(Tv)
        # Infeasible
        @debug "Row $i is primal infeasible"
        ps.status = Trm_PrimalInfeasible
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(Tv)
        ps.solution.y_lower .= zero(Tv)
        ps.solution.y_upper .= zero(Tv)
        ps.solution.s_lower .= zero(Tv)
        ps.solution.s_upper .= zero(Tv)

        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = Sln_Unknown
        ps.solution.dual_status = Sln_InfeasibilityCertificate
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = Tv(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_upper[i] = one(Tv)
        return
    elseif lb > zero(Tv)
        @debug "Row $i is primal infeasible"
        ps.status = Trm_PrimalInfeasible
        ps.updated = true

        # Resize problem
        compute_index_mapping!(ps)
        resize!(ps.solution, ps.nrow, ps.ncol)
        ps.solution.x .= zero(Tv)
        ps.solution.y_lower .= zero(Tv)
        ps.solution.y_upper .= zero(Tv)
        ps.solution.s_lower .= zero(Tv)
        ps.solution.s_upper .= zero(Tv)
        
        # Farkas ray: y⁺_i = 1 (any > 0 value works)
        ps.solution.primal_status = Sln_Unknown
        ps.solution.dual_status = Sln_InfeasibilityCertificate
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = true
        ps.solution.z_primal = ps.solution.z_dual = Tv(Inf)
        i_ = ps.new_con_idx[i]
        ps.solution.y_lower[i_] = one(Tv)
        return
    else
        push!(ps.ops, EmptyRow(i, zero(Tv)))
    end

    # Book-keeping
    ps.updated = true
    ps.rowflag[i] = false
    ps.nrow -= 1
end

function postsolve!(sol::Solution{Tv}, op::EmptyRow{Tv}) where{Tv}
    sol.y_lower[op.i] = pos_part(op.y)
    sol.y_upper[op.i] = neg_part(op.y)
    return nothing
end