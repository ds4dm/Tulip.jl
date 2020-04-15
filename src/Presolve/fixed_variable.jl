struct FixedVariable{Tv} <: PresolveTransformation{Tv}
    j::Int  # variable index
    x::Tv  # primal value
    c::Tv  # current objective coeff
    col::Col{Tv}  # current column
end

function remove_fixed_variable!(lp::PresolveData, j::Int)
    lp.colflag[j] || return nothing  # Column was already removed
    
    # Check bounds
    lb, ub = lp.lcol[j], lp.ucol[j]

    # TODO: use tolerance
    if lb == ub
        @debug "Fix variable $j to $lb"

        col = lp.pb0.acols[j]
        cj = lp.obj[j]

        # Remove column
        lp.colflag[j] = false
        lp.ncol -= 1
        lp.updated = true

        # TODO: make this more efficient
        push!(lp.ops, FixedVariable(j, lb, cj, Col(
            [i for i in col.nzind if lp.rowflag[i]],
            [aij for (i, aij) in zip(col.nzind, col.nzval) if lp.rowflag[i]]
        )))

        # Update objective constant
        lp.obj0 += cj * lb

        # Update rows
        for (i, aij) in zip(col.nzind, col.nzval)
            lp.rowflag[i] || continue  # This row is no longer in the problem

            # Update row bounds
            lp.lrow[i] -= aij * lb
            lp.urow[i] -= aij * lb

            lp.nzrow[i] -= 1

            # Check row
            if lp.nzrow[i] == 0
                remove_empty_row!(lp, i)
            elseif lp.nzrow == 1
                push!(lp.row_singletons, i)
            end
        end  # row update
    end

    # Done
    return nothing
end

function postsolve!(sol::Solution{Tv}, op::FixedVariable{Tv}) where{Tv}
    sol.x[op.j] = op.x
    s = sol.is_dual_ray ? zero(Tv) : op.c
    for (i, aij) in zip(op.col.nzind, op.col.nzval)
        s -= aij * (sol.y_lower[i] - sol.y_upper[i])
    end
    sol.s_lower[op.j] = pos_part(s)
    sol.s_upper[op.j] = neg_part(s)
    return nothing
end