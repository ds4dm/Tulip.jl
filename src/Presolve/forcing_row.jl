struct ForcingRow{T} <: PresolveTransformation{T}
    i::Int  # Row index
    at_lower::Bool  # Whether row is forced to its lower bound (false means upper)
    row::Row{T}  # Row
    cols::Vector{Col{T}}  # Columns of variables in forcing row
    xs::Vector{T}  # Primal values of variables in the row (upper or lower bound)
    cs::Vector{T}  # Objective coeffs of variables in forcing row
end

struct DominatedRow{T} <: PresolveTransformation{T}
    i::Int  # Row index
end

function remove_forcing_row!(ps::PresolveData{T}, i::Int) where{T}
    ps.rowflag[i] || return
    ps.nzrow[i] == 1 && return  # skip row singletons 

    # Implied row bounds
    row = ps.pb0.arows[i]
    l_ = u_ = zero(T)
    for (j, aij) in zip(row.nzind, row.nzval)
        ps.colflag[j] || continue
        if signbit(aij)
            l_ += aij * ps.ucol[j]
            u_ += aij * ps.lcol[j]
        else
            l_ += aij * ps.lcol[j]
            u_ += aij * ps.ucol[j]
        end

        isfinite(l_) || isfinite(u_) || break  # infinite bounds
    end

    l, u = ps.lrow[i], ps.urow[i]
    if l <= l_ <= u_ <= u
        # Constraint is dominated
        @debug "Row $i is dominated"
        ps.rowflag[i] = false
        ps.updated = true
        ps.nrow -= 1

        push!(ps.ops, DominatedRow{T}(i))
        # Update column non-zeros
        for (j, aij) in zip(row.nzind, row.nzval)
            ps.colflag[j] || continue
            ps.nzcol[j] -= !iszero(aij)
        end

    elseif l_ == u
        # Row is forced to its upper bound
        @debug "Row $i is forced to its upper bound"  (l, u) (l_, u_)
        # Record tranformation
        row_ = Row(
            [  j for (j, aij) in zip(row.nzind, row.nzval) if ps.colflag[j]],
            [aij for (j, aij) in zip(row.nzind, row.nzval) if ps.colflag[j]]
        )
        cols_ = Col{T}[]
        xs = T[]
        cs = T[]
        for (j, aij) in zip(row.nzind, row.nzval)
            ps.colflag[j] || continue
            # Extract column j
            col = ps.pb0.acols[j]

            col_ = Col{T}(Int[], T[])

            # Fix variable to its bound
            # TODO: put this in a function and mutualize with fixed variables
            if aij > 0
                # Set xj to its lower bound
                xj_ = ps.lcol[j]
            else
                # Set xj to its upper bound
                xj_ = ps.ucol[j]
            end

            for (k, akj) in zip(col.nzind, col.nzval)
                ps.rowflag[k] || continue

                # Update column j
                push!(col_.nzind, k)
                push!(col_.nzval, akj)

                # Update row k
                ps.nzrow[k] -= 1
                ps.lrow[k] -= akj * xj_
                ps.urow[k] -= akj * xj_

                # ps.nzrow[k] == 0 && remove_empty_row!(ps, k)
                ps.nzrow[k] == 1 && push!(ps.row_singletons, k)
            end
            
            cj = ps.obj[j]
            push!(cols_, col_)
            push!(xs, xj_)
            push!(cs, cj)

            # Remove variable from problem
            ps.colflag[j] = false
            ps.ncol -= 1
        end

        # Record transformation
        push!(ps.ops, ForcingRow(i, true, row_, cols_, xs, cs))

        # Book-keeping
        ps.rowflag[i] = false
        ps.nrow -= 1
        ps.updated = true

    elseif u_ == l
        # Row is forced to its lower bound
        @debug "Row $i is forced to its lower bound" (l, u) (l_, u_)
        # Record tranformation
        row_ = Row(
            [  j for (j, aij) in zip(row.nzind, row.nzval) if ps.colflag[j]],
            [aij for (j, aij) in zip(row.nzind, row.nzval) if ps.colflag[j]]
        )
        cols_ = Col{T}[]
        xs = T[]
        cs = T[]
        for (j, aij) in zip(row.nzind, row.nzval)
            ps.colflag[j] || continue
            # Extract column j
            col = ps.pb0.acols[j]

            col_ = Col{T}(Int[], T[])

            # Fix variable to its bound
            # TODO: put this in a function and mutualize with fixed variables
            if aij > 0
                # Set xj to its upper bound
                xj_ = ps.ucol[j]
            else
                # Set xj to its lower bound
                xj_ = ps.lcol[j]
            end

            for (k, akj) in zip(col.nzind, col.nzval)
                ps.rowflag[k] || continue

                # Update column j
                push!(col_.nzind, k)
                push!(col_.nzval, akj)

                # Update row k
                ps.nzrow[k] -= 1
                ps.lrow[k] -= akj * xj_
                ps.urow[k] -= akj * xj_

                # ps.nzrow[k] == 0 && remove_empty_row!(ps, k)
                ps.nzrow[k] == 1 && push!(ps.row_singletons, k)
            end
            
            cj = ps.obj[j]
            push!(cols_, col_)
            push!(xs, xj_)
            push!(cs, cj)

            # Remove variable from problem
            ps.colflag[j] = false
            ps.ncol -= 1
        end

        # Record transformation
        push!(ps.ops, ForcingRow(i, false, row_, cols_, xs, cs))

        # Book-keeping
        ps.rowflag[i] = false
        ps.nrow -= 1
        ps.updated = true
    end
    # TODO: handle infeasible row cases

    return nothing
end

function postsolve!(sol::Solution{T}, op::DominatedRow{T}) where{T}
    sol.y_lower[op.i] = zero(T)
    sol.y_upper[op.i] = zero(T)
    return nothing
end

# TODO: postsolve of forcing rows
function postsolve!(sol::Solution{T}, op::ForcingRow{T}) where{T}

    # Primal
    for (j, xj) in zip(op.row.nzind, op.xs)
        sol.x[j] = xj
    end

    # Dual
    z = similar(op.cs)
    for (l, (j, cj, col)) in enumerate(zip(op.row.nzind, op.cs, op.cols))
        z[l] = cj
        for (k, akj) in zip(col.nzind, col.nzval)
            z[l] -= akj * (sol.y_lower[k] - sol.y_upper[k])
        end
    end

    # First, compute yi
    y = op.at_lower ? maximum(z ./ op.row.nzval) : minimum(z ./ op.row.nzval)
    sol.y_lower[op.i] = pos_part(y)
    sol.y_upper[op.i] = neg_part(y)

    # Extract reduced costs
    for (j, aij, zj) in zip(op.row.nzind, op.row.nzval, z)
        s = zj - aij * y
        sol.s_lower[j] = pos_part(s)
        sol.s_upper[j] = neg_part(s)
    end
    return nothing
end
