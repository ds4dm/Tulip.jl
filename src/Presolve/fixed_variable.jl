struct FixedVariable{T} <: PresolveTransformation{T}
    j::Int  # variable index
    x::T  # primal value
    c::T  # current objective coeff
    col::Col{T}  # current column
end

function remove_fixed_variable!(ps::PresolveData, j::Int)
    ps.colflag[j] || return nothing  # Column was already removed
    
    # Check bounds
    lb, ub = ps.lcol[j], ps.ucol[j]

    # TODO: use tolerance
    if lb == ub
        @debug "Fix variable $j to $lb"

        col = ps.pb0.acols[j]
        cj = ps.obj[j]

        # Remove column
        ps.colflag[j] = false
        ps.ncol -= 1
        ps.updated = true

        # TODO: make this more efficient
        push!(ps.ops, FixedVariable(j, lb, cj, Col(
            [i for i in col.nzind if ps.rowflag[i]],
            [aij for (i, aij) in zip(col.nzind, col.nzval) if ps.rowflag[i]]
        )))

        # Update objective constant
        ps.obj0 += cj * lb

        # Update rows
        for (i, aij) in zip(col.nzind, col.nzval)
            ps.rowflag[i] || continue  # This row is no longer in the problem
            iszero(aij) && continue  # Skip if coefficient is zero

            # Update row bounds
            ps.lrow[i] -= aij * lb
            ps.urow[i] -= aij * lb

            ps.nzrow[i] -= 1

            # Check row
            if ps.nzrow[i] == 0
                remove_empty_row!(ps, i)
            elseif ps.nzrow == 1
                push!(ps.row_singletons, i)
            end
        end  # row update
    end

    # Done
    return nothing
end

function postsolve!(sol::Solution{T}, op::FixedVariable{T}) where{T}
    sol.x[op.j] = op.x
    s = sol.is_dual_ray ? zero(T) : op.c
    for (i, aij) in zip(op.col.nzind, op.col.nzval)
        s -= aij * (sol.y_lower[i] - sol.y_upper[i])
    end
    sol.s_lower[op.j] = pos_part(s)
    sol.s_upper[op.j] = neg_part(s)
    return nothing
end