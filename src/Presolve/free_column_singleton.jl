struct FreeColumnSingleton{Tv} <: PresolveTransformation{Tv}
    i::Int  # Row index
    j::Int  # Column index
    l::Tv  # Row lower bound
    u::Tv  # Row upper bound
    aij::Tv
    y::Tv  # Dual variable
    row::Row{Tv}
end

function remove_free_column_singleton!(ps::PresolveData{Tv}, j::Int) where{Tv}

    ps.colflag[j] && ps.nzcol[j] == 1 || return nothing  # only column singletons

    col = ps.pb0.acols[j]

    # Find non-zero index
    # TODO: put this in a function
    nz = 0
    i, aij = 0, zero(Tv)
    for (i_, a_) in zip(col.nzind, col.nzval)
        if ps.rowflag[i_]
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

    row = ps.pb0.arows[i]
    lr, ur = ps.lrow[i], ps.urow[i]

    # Detect if xj is implied free
    l, u = ps.lcol[j], ps.ucol[j]
    if isfinite(l) || isfinite(u)
        # Not a free variable, compute implied bounds
        if aij > zero(Tv)
            l_, u_ = lr, ur
            for (k, aik) in zip(row.nzind, row.nzval)
                (ps.colflag[k] && k != j) || continue
                # Update bounds
                if aik > 0
                    l_ -= aik * ps.ucol[k]
                    u_ -= aik * ps.lcol[k]
                else
                    l_ -= aik * ps.lcol[k]
                    u_ -= aik * ps.ucol[k]
                end
            end
            l_ /= aij
            u_ /= aij
        else
            l_, u_ = ur, lr
            for (k, aik) in zip(row.nzind, row.nzval)
                (ps.colflag[k] && k != j) || continue
                # Update bounds
                if aik > 0
                    l_ -= aik * ps.lcol[k]
                    u_ -= aik * ps.ucol[k]
                else
                    l_ -= aik * ps.ucol[k]
                    u_ -= aik * ps.lcol[k]
                end
            end
            l_ /= aij
            u_ /= aij
        end
        @debug """Column singleton $j
            Original bounds: [$l, $u]
             Implied bounds: [$(l_), $(u_)]
        """
        l <= l_ <= u_ <= u || return nothing  # Not implied free
    end

    # Remove (implied) free column
    @debug "Removing (implied) free column singleton $j"
    y = ps.obj[j] / aij  # dual of row i

    # Update objective
    ps.obj0 += (y >= zero(Tv)) ? y * lr : y * ur
    row_ = Row{Tv}(Int[], Tv[])
    for (j_, aij_) in zip(row.nzind, row.nzval)
        ps.colflag[j_] && (j_ != j) || continue

        push!(row_.nzind, j_)
        push!(row_.nzval, aij_)
        ps.obj[j_] -= y * aij_

        # Update number of non-zeros in column
        ps.nzcol[j_] -= 1
    end

    # Book-keeping
    push!(ps.ops, FreeColumnSingleton(i, j, lr, ur, aij, y, row_))
    ps.rowflag[i] = false  # remove row
    ps.colflag[j] = false  # remove column
    ps.nrow -= 1
    ps.ncol -= 1
    ps.updated = true

    return nothing
end

function postsolve!(sol::Solution{Tv}, op::FreeColumnSingleton{Tv}) where{Tv}
    # Dual
    y = op.y
    sol.y_lower[op.i] = pos_part(y)
    sol.y_upper[op.i] = neg_part(y)
    sol.s_lower[op.j] = zero(Tv)
    sol.s_upper[op.j] = zero(Tv)

    # Primal
    sol.x[op.j] = sol.is_primal_ray ? zero(Tv) : (y >= zero(Tv) ? op.l : op.u)
    for (k, aik) in zip(op.row.nzind, op.row.nzval)
        sol.x[op.j] -= aik * sol.x[k]
    end
    sol.x[op.j] /= op.aij
    
    return nothing
end