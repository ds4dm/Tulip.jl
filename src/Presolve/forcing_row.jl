struct ForcingRow{Tv} <: PresolveTransformation{Tv}
    i::Int  # Row index
    at_lower::Bool  # Whether row is forced to its lower bound (false means upper)
    row::Row{Tv}  # Row
    cols::Vector{Col{Tv}}  # Columns of variables in forcing row
    xs::Vector{Tv}  # Primal values of variables in the row (upper or lower bound)
    cs::Vector{Tv}  # Objective coeffs of variables in forcing row
end

struct DominatedRow{Tv} <: PresolveTransformation{Tv}
    i::Int  # Row index
end

function remove_forcing_row!(lp::PresolveData{Tv}, i::Int) where{Tv}
    lp.rowflag[i] || return
    lp.nzrow[i] == 1 && return  # skip row singletons 

    # Implied row bounds
    row = lp.pb0.arows[i]
    l_ = u_ = zero(Tv)
    for (j, aij) in zip(row.nzind, row.nzval)
        lp.colflag[j] || continue
        if aij < zero(Tv)
            l_ += aij * lp.ucol[j]
            u_ += aij * lp.lcol[j]
        else
            l_ += aij * lp.lcol[j]
            u_ += aij * lp.ucol[j]
        end

        isfinite(l_) || isfinite(u_) || break  # infinite bounds
    end

    l, u = lp.lrow[i], lp.urow[i]
    if l <= l_ <= u_ <= u
        # Constraint is dominated
        @debug "Row $i is dominated"
        lp.rowflag[i] = false
        lp.updated = true
        lp.nrow -= 1

        push!(lp.ops, DominatedRow{Tv}(i))
        # Update column non-zeros
        for (j, aij) in zip(row.nzind, row.nzval)
            lp.colflag[j] || continue
            lp.nzcol[j] -= !iszero(aij)
        end

    elseif l_ == u
        # Row is forced to its upper bound
        @debug "Row $i is forced to its upper bound"  (l, u) (l_, u_)
        # Record tranformation
        row_ = Row(
            [  j for (j, aij) in zip(row.nzind, row.nzval) if lp.colflag[j]],
            [aij for (j, aij) in zip(row.nzind, row.nzval) if lp.colflag[j]]
        )
        cols_ = Col{Tv}[]
        xs = Tv[]
        cs = Tv[]
        for (j, aij) in zip(row.nzind, row.nzval)
            lp.colflag[j] || continue
            # Extract column j
            col = lp.pb0.acols[j]

            col_ = Col{Tv}(Int[], Tv[])

            # Fix variable to its bound
            # TODO: put this in a function and mutualize with fixed variables
            if aij > 0
                # Set xj to its lower bound
                xj_ = lp.lcol[j]
            else
                # Set xj to its upper bound
                xj_ = lp.ucol[j]
            end

            for (k, akj) in zip(col.nzind, col.nzval)
                lp.rowflag[k] || continue

                # Update column j
                push!(col_.nzind, k)
                push!(col_.nzval, akj)

                # Update row k
                lp.nzrow[k] -= 1
                lp.lrow[k] -= akj * xj_
                lp.urow[k] -= akj * xj_

                # lp.nzrow[k] == 0 && remove_empty_row!(lp, k)
                lp.nzrow[k] == 1 && push!(lp.row_singletons, k)
            end
            
            cj = lp.obj[j]
            push!(cols_, col_)
            push!(xs, xj_)
            push!(cs, cj)

            # Remove variable from problem
            lp.colflag[j] = false
            lp.ncol -= 1
        end

        # Record transformation
        push!(lp.ops, ForcingRow(i, true, row_, cols_, xs, cs))

        # Book-keeping
        lp.rowflag[i] = false
        lp.nrow -= 1
        lp.updated = true

    elseif u_ == l
        # Row is forced to its lower bound
        @debug "Row $i is forced to its lower bound" (l, u) (l_, u_)
        # Record tranformation
        row_ = Row(
            [  j for (j, aij) in zip(row.nzind, row.nzval) if lp.colflag[j]],
            [aij for (j, aij) in zip(row.nzind, row.nzval) if lp.colflag[j]]
        )
        cols_ = Col{Tv}[]
        xs = Tv[]
        cs = Tv[]
        for (j, aij) in zip(row.nzind, row.nzval)
            lp.colflag[j] || continue
            # Extract column j
            col = lp.pb0.acols[j]

            col_ = Col{Tv}(Int[], Tv[])

            # Fix variable to its bound
            # TODO: put this in a function and mutualize with fixed variables
            if aij > 0
                # Set xj to its upper bound
                xj_ = lp.ucol[j]
            else
                # Set xj to its lower bound
                xj_ = lp.lcol[j]
            end

            for (k, akj) in zip(col.nzind, col.nzval)
                lp.rowflag[k] || continue

                # Update column j
                push!(col_.nzind, k)
                push!(col_.nzval, akj)

                # Update row k
                lp.nzrow[k] -= 1
                lp.lrow[k] -= akj * xj_
                lp.urow[k] -= akj * xj_

                # lp.nzrow[k] == 0 && remove_empty_row!(lp, k)
                lp.nzrow[k] == 1 && push!(lp.row_singletons, k)
            end
            
            cj = lp.obj[j]
            push!(cols_, col_)
            push!(xs, xj_)
            push!(cs, cj)

            # Remove variable from problem
            lp.colflag[j] = false
            lp.ncol -= 1
        end

        # Record transformation
        push!(lp.ops, ForcingRow(i, false, row_, cols_, xs, cs))

        # Book-keeping
        lp.rowflag[i] = false
        lp.nrow -= 1
        lp.updated = true
    end
    # TODO: handle infeasible row cases

    return nothing
end

function postsolve!(sol::Solution{Tv}, op::DominatedRow{Tv}) where{Tv}
    sol.y_lower[op.i] = zero(Tv)
    sol.y_upper[op.i] = zero(Tv)
    return nothing
end

# TODO: postsolve of forcing rows
function postsolve!(sol::Solution{Tv}, op::ForcingRow{Tv}) where{Tv}

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