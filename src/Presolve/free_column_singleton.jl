struct FreeColumnSingleton{Tv} <: PresolveTransformation{Tv}
    i::Int  # Row index
    j::Int  # Column index
    l::Tv  # Row lower bound
    u::Tv  # Row upper bound
    aij::Tv
    y::Tv  # Dual variable
    row::Row{Tv}
end

function remove_free_column_singleton!(lp::PresolveData{Tv}, j::Int) where{Tv}

    lp.colflag[j] && lp.nzcol[j] == 1 || return nothing  # only column singletons

    col = lp.pb0.acols[j]

    # Find non-zero index
    nz = 0
    i, aij = 0, zero(Tv)
    for (i_, a_) in zip(col.nzind, col.nzval)
        if lp.rowflag[i_]
            nz += 1; nz <= 1 || break

            i = i_
            aij = a_
        end
    end

    # Not a singleton
    nz == 1 || (@error "Expected singletons but column $j has $nz non-zeros"; return nothing)

    if iszero(aij) || iszero(i)
        return nothing  # column is actually empty
    end

    row = lp.pb0.arows[i]
    lr, ur = lp.lrow[i], lp.urow[i]

    # Detect if xj is implied free
    l, u = lp.lcol[j], lp.ucol[j]
    if isfinite(l) || isfinite(u)
        # Not a free variable, compute implied bounds
        if aij > zero(Tv)
            l_, u_ = lr, ur
            for (k, aik) in zip(row.nzind, row.nzval)
                lp.colflag[k] || continue

                # Update bounds
                if aik > 0
                    l_ -= aik * lp.ucol[k]
                    u_ -= aik * lp.lcol[k]
                else
                    l_ -= aik * lp.lcol[k]
                    u_ -= aik * lp.ucol[k]
                end
            end
            l_ /= aij
            u_ /= aij
        else
            l_, u_ = ur, lr
            for (k, aik) in zip(row.nzind, row.nzval)
                lp.colflag[k] || continue

                # Update bounds
                if aik > 0
                    l_ -= aik * lp.lcol[k]
                    u_ -= aik * lp.ucol[k]
                else
                    l_ -= aik * lp.ucol[k]
                    u_ -= aik * lp.lcol[k]
                end
            end
            l_ /= aij
            u_ /= aij
        end

        l <= l_ <= u_ <= u || return nothing  # Not implied free
    end

    # Remove (implied) free column
    @debug "Remove free column singleton $j"
    y = lp.obj[j] / aij  # dual of row i

    # Update objective
    lp.obj0 += (y >= zero(Tv)) ? y * lr : y * ur
    row_ = Row{Tv}(Int[], Tv[])
    for (j_, aij_) in zip(row.nzind, row.nzval)
        lp.colflag[j_] && (j_ != j) || continue

        push!(row_.nzind, j_)
        push!(row_.nzval, aij_)
        lp.obj[j_] -= y * aij_

        # Update number of non-zeros in column
        lp.nzcol[j_] -= 1
    end

    # Book-keeping
    push!(lp.ops, FreeColumnSingleton(i, j, lr, ur, aij, y, row_))
    lp.rowflag[i] = false  # remove row
    lp.colflag[j] = false  # remove column
    lp.nrow -= 1
    lp.ncol -= 1
    lp.updated = true

    return nothing
end

function postsolve!(sol::Solution{Tv}, sol_::Solution{Tv}, op::FreeColumnSingleton{Tv}) where{Tv}
    # Dual
    y = op.y
    sol.y_lower[op.i] = pos_part(y)
    sol.y_upper[op.i] = neg_part(y)
    sol.s_lower[op.j] = zero(Tv)
    sol.s_upper[op.j] = zero(Tv)

    # Primal
    sol.x[op.j] = y >= zero(Tv) ? op.l : op.u
    for (k, aik) in zip(op.row.nzind, op.row.nzval)
        sol.x[op.j] -= aik * sol.x[k]
    end
    sol.x[op.j] /= op.aij
    
    return nothing
end