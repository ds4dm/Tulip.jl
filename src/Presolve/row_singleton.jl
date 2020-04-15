struct RowSingleton{Tv} <: PresolveTransformation{Tv}
    i::Int  # Row index
    j::Int  # Column index
    aij::Tv  # Row coefficient
    force_lower::Bool  # Whether row was forcing the lower bound
    force_upper::Bool  # Whether row was forcing the upper bound
end


function remove_row_singleton!(lp::PresolveData{Tv}, i::Int) where{Tv}
    # Sanity checks
    (lp.rowflag[i] && lp.nzrow[i] == 1) || return nothing

    row = lp.pb0.arows[i]

    # Find non-zero coefficient and its column index
    nz = 0
    j, aij = 0, zero(Tv)
    for (j_, aij_) in zip(row.nzind, row.nzval)
        if lp.colflag[j_] && !iszero(aij_)
            nz += 1; nz <= 1 || break  # not a row singleton

            j = j_
            aij = aij_
        end
    end

    if nz > 1
        @error "Row $i was marked as row singleton but has $nz non-zeros"
        return nothing
    end

    if iszero(aij)
        # Row is actually empty
        # It will be removed at the next forcing constraints check
        return nothing
    end

    # Compute implied bounds
    if aij > zero(Tv)
        l = lp.lrow[i] / aij
        u = lp.urow[i] / aij
    else
        l = lp.urow[i] / aij
        u = lp.lrow[i] / aij
    end

    # Compare to variable bounds
    # TODO: what if bounds are incompatible?
    lb, ub = lp.lcol[j], lp.ucol[j]
    force_lower = (l >= lb)
    force_upper = (u <= ub)
    if l >= lb
        lp.lcol[j] = l
    end
    if u <= ub
        lp.ucol[j] = u
    end

    # Book-keeping
    push!(lp.ops, RowSingleton(i, j, aij, force_lower, force_upper))
    lp.rowflag[i] = false
    lp.updated = true
    lp.nrow -= 1
    lp.nzcol[j] -= 1

    # Check if we created a fixed/empty column
    if lp.lcol[j] == lp.ucol[j]
        remove_fixed_variable!(lp, j)
        return nothing
    end
    if lp.nzcol[j] == 0
        # TODO: remove empty column
    elseif lp.nzcol[j] == 1
        # TODO: track column singleton
    end

    return nothing
end

function postsolve!(sol::Solution{Tv}, op::RowSingleton{Tv}) where{Tv}

    if op.force_lower
        if op.aij > zero(Tv)
            sol.y_lower[op.i] = sol.s_lower[op.j] / op.aij
        else
            sol.y_upper[op.i] = sol.s_lower[op.j] / abs(op.aij)
        end
        sol.s_lower[op.j] = zero(Tv)
    end
    if op.force_upper
        if op.aij > zero(Tv)
            sol.y_upper[op.i] = sol.s_upper[op.j] / op.aij
        else
            sol.y_lower[op.i] = sol.s_upper[op.j] / abs(op.aij)
        end
        sol.s_upper[op.j] = zero(Tv)
    end

    return nothing
end