struct RowSingleton{T} <: PresolveTransformation{T}
    i::Int  # Row index
    j::Int  # Column index
    aij::T  # Row coefficient
    force_lower::Bool  # Whether row was forcing the lower bound
    force_upper::Bool  # Whether row was forcing the upper bound
end


function remove_row_singleton!(ps::PresolveData{T}, i::Int) where{T}
    # Sanity checks
    (ps.rowflag[i] && ps.nzrow[i] == 1) || return nothing

    row = ps.pb0.arows[i]

    # Find non-zero coefficient and its column index
    nz = 0
    j, aij = 0, zero(T)
    for (j_, aij_) in zip(row.nzind, row.nzval)
        if ps.colflag[j_] && !iszero(aij_)
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
    if !signbit(aij)
        l = ps.lrow[i] / aij
        u = ps.urow[i] / aij
    else
        l = ps.urow[i] / aij
        u = ps.lrow[i] / aij
    end

    # Compare to variable bounds
    # TODO: what if bounds are incompatible?
    lb, ub = ps.lcol[j], ps.ucol[j]
    force_lower = (l >= lb)
    force_upper = (u <= ub)
    if l >= lb
        ps.lcol[j] = l
    end
    if u <= ub
        ps.ucol[j] = u
    end

    # Book-keeping
    push!(ps.ops, RowSingleton(i, j, aij, force_lower, force_upper))
    ps.rowflag[i] = false
    ps.updated = true
    ps.nrow -= 1
    ps.nzcol[j] -= 1

    # Check if we created a fixed/empty column
    if ps.lcol[j] == ps.ucol[j]
        remove_fixed_variable!(ps, j)
        return nothing
    end
    if ps.nzcol[j] == 0
        # TODO: remove empty column
    elseif ps.nzcol[j] == 1
        # TODO: track column singleton
    end

    return nothing
end

function postsolve!(sol::Solution{T}, op::RowSingleton{T}) where{T}

    if op.force_lower
        if op.aij > zero(T)
            sol.y_lower[op.i] = sol.s_lower[op.j] / op.aij
        else
            sol.y_upper[op.i] = sol.s_lower[op.j] / abs(op.aij)
        end
        sol.s_lower[op.j] = zero(T)
    end
    if op.force_upper
        if op.aij > zero(T)
            sol.y_upper[op.i] = sol.s_upper[op.j] / op.aij
        else
            sol.y_lower[op.i] = sol.s_upper[op.j] / abs(op.aij)
        end
        sol.s_upper[op.j] = zero(T)
    end

    return nothing
end
