"""
    readmps(filename)

    Parse an MPS file, and return a linear program in standard form.
"""
function readmps!(m::Model{Tv}, fname::String) where{Tv<:Real}

    # First, empty model
    empty!(m)

    # Now, parse MPS file
    section = ""
    d = Dict{String,Any}()

    con2idx = Dict{String, ConstrId}()  # name <-> index correspondence for rows
    var2idx = Dict{String, VarId}()     # name <-> index correspondence for columns

    lb = Dict{VarId, Tv}()  # Lower bounds on variables
    ub = Dict{VarId, Tv}()  # Upper bounds on variables

    d["nrows"] = 0
    d["ncols"] = 0
    d["rhs"] = ""
    d["bound"] = ""
    d["range"] = ""

    open(fname) do f
        for ln in eachline(f)

            fields = parseline(ln)

            # pass empty lines or comments
            if length(fields) == 0
                continue
            end

            # check for indicators
            if ln[1] != ' '
                if fields[1] == "NAME"
                    # extract problem name
                    if length(fields) >= 2
                        # second field is name, rest is ignored
                        m.name = fields[2]
                    else
                        #  empty name
                        m.name = ""
                    end
                    continue
                end
                # update section value and pass to next line
                section = fields[1]
                continue
            end

            if section == "ROWS"
                parsemps_rows!(m, ln, fields, con2idx)

            elseif section == "COLUMNS"
                parsemps_columns!(m, ln, fields, var2idx, con2idx, lb, ub)

            elseif section == "RHS"
                parsemps_rhs!(m, ln, fields, con2idx, d)

            elseif section == "RANGES"
                parsemps_ranges!(m, ln, fields, d, con2idx)

            elseif section == "BOUNDS"
                parsemps_bounds!(m, ln, fields, d, var2idx, lb, ub)

            elseif section == "ENDATA"
                # end of file
            else
                @warn "UNEXPECTED FORMAT PARSING\t$(ln)"
            end

        end

        # End of file reached
    end

    # Set variable bounds
    for (vidx, var) in m.pbdata_raw.vars
        set_bounds!(var, lb[vidx], ub[vidx])
    end

    return m
end


"""
    parseline(ln)
"""
function parseline(ln)
    s = []
    if length(ln) == 0 || ln[1] == '*'
        # empty line or comment line
        return s
    elseif ln[1] != ' '
        # indicator row
        s = split(ln)
        return s
    end

    # fill-in with spaces to get a 80-long line
    if length(ln) < 80
        ln_ = ln * (" "^(80-length(ln)))
    end

    fields = [
        ln_[2:3],
        ln_[5:12],
        ln_[15:22],
        ln_[25:36],
        ln_[40:47],
        ln_[50:61]
    ]

    # parse each field from right to left
    # empty fields on the right are ignored
    for f in fields[end:-1:1]
        # remove left-trailing and right-trailing spaces
        f_ = strip(f)  # f_ may be empty

        if length(f_) > 0 || length(s) > 0
            pushfirst!(s, f_)
        end
    end

    return s
end


function parsemps_rows!(m::Model{Tv}, ln, fields, con2idx) where{Tv<:Real}
    # parse a ROWS line of the MPS file
    # `ln` is a string that contains the current line
    # `fields` is an array such that fields == split(ln)
    
    # current line should contain exactly 2 fields
    if length(fields) != 2
        # input error, current line is ignored
        @warn "INPUT ERROR:\t$(ln)"
        return nothing
    end

    # parse line
    # First field can be either of:
    #   "N" -> row corresponds to objective
    #   "E" -> equality constraint
    #   "L" -> less-or-equal inequality constraint
    #   "G" -> greater-or-equal inequality constraint
    if fields[1] == "N"
        # objective
        con2idx[fields[2]] = ConstrId(0)
    elseif fields[1] == "E" || fields[1] == "L" || fields[1] == "G"
        # constraint
        if fields[1] == "E"
            lb, ub = zero(Tv), zero(Tv)
        elseif fields[1] == "L"
            lb, ub = Tv(-Inf), zero(Tv)
        elseif fields[1] == "G"
            lb, ub = zero(Tv), Tv(Inf)
        end

        # add constraint to model
        ridx = add_constraint!(m, String(fields[2]), lb, ub, VarId[], Float64[])
        con2idx[String(fields[2])] = ridx
    else
        # input error, current line is ignored
        @warn "INPUT ERROR:\t$(ln)"
        return nothing
    end
end


function parsemps_columns!(m::Model{Tv}, ln, fields, var2idx, con2idx, lb, ub) where{Tv<:Real}
    # parse a ROWS line of the MPS file
    # `ln` is a string that contains the current line
    # `fields` contains the 6 fields of that line
    
    # current line should contain at least three fields
    if length(fields) < 4
        @warn "COLUMN LINE TOO SHORT|$(ln)"
        return nothing
    end

    # First field is empty, second field is variable's name
    cname = fields[2]
    if !haskey(var2idx, cname)
        # d["ncols"] += 1
        # var2idx[cname] = d["ncols"]  # Add column name to the pool
        # Create new variable
        vidx = add_variable!(m, String(cname), zero(Tv), zero(Tv), Tv(Inf))
        var2idx[cname] = vidx
        lb[vidx] = zero(Tv)
        ub[vidx] = Tv(Inf)
    end

    # Second and third fields are
    # the row's name, and the corresponding coefficient
    rname = fields[3]  # row name
    if !haskey(con2idx, rname)
        # current row not in list of rows, current line ignored
        @warn "COL ERROR UNKNOWN ROW $(rname)|$(ln)"
        return nothing
    end
    coeff = parse(Float64, fields[4])  # coefficient value
    if con2idx[rname].uuid == 0
        # objective coefficient
        # push!(obj_col, var2idx[cname])
        # push!(obj_val, coeff)
        set_obj_coeff!(m.pbdata_raw.vars[var2idx[cname]], coeff)
    else
        # constraint coefficient
        set_coeff!(m.pbdata_raw, var2idx[cname], con2idx[rname], coeff)
        # push!(coeffs_row, con2idx[rname])
        # push!(coeffs_val, coeff)
    end

    # optional other fields
    if length(fields) >= 6
        rname = fields[5]  # row name
        if !haskey(con2idx, rname)
            # current row not in list of rows, current line ignored
            @warn "COL ERROR UNKNOWN ROW $(rname)|$(ln)"
            return nothing
        end
        coeff = parse(Float64, fields[6])  # coefficient value
        if con2idx[rname].uuid == 0
            # objective coefficient
            # push!(obj_col, var2idx[cname])
            # push!(obj_val, coeff)
            set_obj_coeff!(m.pbdata_raw.vars[var2idx[cname]], coeff)
        else
            # constraint coefficient
            # push!(coeffs_col, var2idx[cname])
            # push!(coeffs_row, con2idx[rname])
            # push!(coeffs_val, coeff)
            set_coeff!(m.pbdata_raw, var2idx[cname], con2idx[rname], coeff)
        end
    end

    return nothing
end


function parsemps_rhs!(m::Model{Tv}, ln, fields, con2idx, d) where{Tv<:Real}
    # parse a line of the RHS section
    # `ln` is the current line
    # `fields` contains the fields of that line


    # Check for too short lines
    if length(fields) < 4
        # input error, current line is ignored
        @warn "RHS LINE TOO SHORT|$(ln)"
        return nothing
    end

    # First field is empty
    # Second field is RHS name (may be empty)
    rhsname = fields[2]
    if d["rhs"] == ""
        # first time RHS is read
        d["rhs"] = rhsname
    elseif d["rhs"] != rhsname
        # other RHS, current line is ignored
        return nothing
    end

    
    # parse line
    rname = fields[3]
    if !haskey(con2idx, rname)
        # current row not in list of rows, current line ignored
        @warn "RHS ERROR UNKNOWN ROW $(rname)|$(ln)"
        return nothing
    end
    rval = parse(Float64, fields[4])
    # update index and value
    if con2idx[rname].uuid == 0
        d["obj_offset"] = -rval
    else
        # Check type of constraint and update coefficient accordingly
        (bt, lb, ub) = get_bounds(m.pbdata_raw.constrs[con2idx[rname]])
        if bt == TLP_UP
            # a'x <= b
            set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], Tv(-Inf), rval)
        elseif bt == TLP_LO
            # a'x >= b
            set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], rval, Tv(Inf))
        elseif bt == TLP_FR
            # This should not happen
            error("Got right-hand side for free constraint")
        elseif bt == TLP_FX
            # a'x = b
            set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], rval, rval)
        elseif bt == TLP_RG
            # This should not happen
            error("Got single right-hand side for ranged constraint.")
        end
        # push!(rhs_row, con2idx[rname])
        # push!(rhs_val, rval)
    end

    # optional fields
    if length(fields) >= 6
        rname = fields[5]  # row name
        if !haskey(con2idx, rname)
            # current row not in list of rows, current line ignored
            @warn "RHS ERROR UNKNOWN ROW $(rname)|$(ln)"
            return nothing
        end
        rval = parse(Float64, fields[6])  # coefficient value
        # update index and value
        if con2idx[rname].uuid == 0
            d["obj_offset"] = -rval
        else
            # Check type of constraint and update coefficient accordingly
            (bt, lb, ub) = get_bounds(m.pbdata_raw.constrs[con2idx[rname]])
            if bt == TLP_UP
                # a'x <= b
                set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], Tv(-Inf), rval)
            elseif bt == TLP_LO
                # a'x >= b
                set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], rval, Tv(Inf))
            elseif bt == TLP_FR
                # This should not happen
                error("Got right-hand side for free constraint")
            elseif bt == TLP_FX
                # a'x = b
                set_bounds!(m.pbdata_raw.constrs[con2idx[rname]], rval, rval)
            elseif bt == TLP_RG
                # This should not happen
                error("Got single right-hand side for ranged constraint.")
            end
        end
    end
end


function parsemps_ranges!(m::Model{Tv}, ln, fields, d, con2idx) where{Tv<:Real}
    if length(fields) < 4
        # input error, current line is ignored
        @warn "RNG LINE TOO SHORT|$(ln)"
        return nothing
    end

    rngname = fields[2]
    if d["range"] == ""
        d["range"] = rngname
    elseif d["range"] != rngname
        # other range vector, ignore
        return nothing
    end

    # parse line
    rname = fields[3]
    rval = parse(Float64, fields[4])
    if !haskey(con2idx, rname)
        # unknown row
        @warn "RNG ERROR UNKNOWN ROW $(rname)|$(ln)"
        return nothing
    end
    cidx = con2idx[rname]
    (bt, lb, ub) = get_bounds(m.pbdata_raw.constrs[cidx])
    if bt == TLP_LO
        # `l <= a'x <= Inf` becomes `l <= a'x <= l + |r|`
        set_bounds!(m.pbdata_raw.constrs[cidx], lb, lb + abs(rval))
    elseif bt == TLP_UP
        # `-Inf <= a'x <= u` becomes `u - |r| <= a'x <= u`
        set_bounds!(m.pbdata_raw.constrs[cidx], ub - abs(rval), ub)
    elseif bt == TLP_FX && rval >= 0.0
        # `a'x = b` becomes `b <= a'x <= b + |r|`
        set_bounds!(m.pbdata_raw.constrs[cidx], lb, lb + rval)
    elseif bt == TLP_FX && rval < 0.0
        # `a'x = b` becomes `b - |r| <= a'x <= b`
        set_bounds!(m.pbdata_raw.constrs[cidx], ub - abs(rval), ub)
    else
        error("Unkown row type for RANGES: $bt.")
    end

    if length(fields) >=6
        rname = fields[5]
        rngval = parse(Float64, fields[6])
        if !haskey(con2idx, rname)
            # unknown row
            @warn "RNG ERROR UNKNOWN ROW $(rname)|$(ln)"
            return nothing
        end
        cidx = con2idx[rname]
        (bt, lb, ub) = get_bounds(m.pbdata_raw.constrs[cidx])
        if bt == TLP_LO
            # `l <= a'x <= Inf` becomes `l <= a'x <= l + |r|`
            set_bounds!(m.pbdata_raw.constrs[cidx], lb, lb + abs(rval))
        elseif bt == TLP_UP
            # `-Inf <= a'x <= u` becomes `u - |r| <= a'x <= u`
            set_bounds!(m.pbdata_raw.constrs[cidx], ub - abs(rval), ub)
        elseif bt == TLP_FX && rval >= 0.0
            # `a'x = b` becomes `b <= a'x <= b + |r|`
            set_bounds!(m.pbdata_raw.constrs[cidx], lb, lb + rval)
        elseif bt == TLP_FX && rval < 0.0
            # `a'x = b` becomes `b - |r| <= a'x <= b`
            set_bounds!(m.pbdata_raw.constrs[cidx], ub - abs(rval), ub)
        else
            error("Unkown row type for RANGES: $bt.")
        end
    end
end


function parsemps_bounds!(m::Model{Tv}, ln, fields, d, var2idx, lb, ub) where{Tv<:Real}
    # parse a line of the BOUNDS section

    # Check
    

    if (length(fields) < 3
        || ((length(fields) < 4) 
            && (fields[1] == "LO" || fields[1] == "UP" || fields[1] == "FX"))
    )
        @warn "BND ERROR LINE TOO SHORT|$(ln)"
        return nothing
    end

    # first field should be bound type
    btype = fields[1]
    bname = fields[2]
    cname = fields[3]
    if !haskey(var2idx, cname)
        @warn "BND ERROR UNKNOWN COL $(cname)|$(ln)"
        return nothing
    end

    if d["bound"] == ""
        d["bound"] = bname
    else
        # Some field was already read; if this one is different, skip
        d["bound"] == bname || return nothing
    end

    vidx = var2idx[cname]

    if length(fields) == 3
        if btype == "FR"
            # -Inf < x < Inf
            lb[vidx] = -Inf
            ub[vidx] = Inf
            return
        elseif btype == "MI"
            # -Inf < x 
            lb[vidx] = -Inf

        elseif btype == "PL"
            # x < Inf
            ub[vidx] = Inf
        
        elseif btype == "BV"
            # x = 0 or 1, x binary
            # Keep bounds but ignore binary requirement
            @warn "Binary variable $cname; recording bounds but ignoring binary."
            lb[vidx] = 0.0
            ub[vidx] = 1.0
        else
            @warn "Unknown bound type: $(btype). Input ignored."
        end
        return
    end

    bval = parse(Float64, fields[4])

    if btype == "LO"
        # b <= x
        lb[vidx] = bval

    elseif btype == "UP"
        # x <= b
        ub[vidx] = bval

    elseif btype == "FX"
        # x = b
        lb[vidx] = bval
        ub[vidx] = bval

    elseif btype == "LI"
        # 0 <= x < Inf, x integer
        # Keep bounds but ignore integer requirement
        @warn "Integer variable $cname; recording bounds but ignoring binary."
        lb[vidx] = bval

    elseif btype == "UI"
        # 0 <= x <= b, x integer
        # Keep bounds but ignore integer requirement
        @warn "Integer variable $cname; recording bounds but ignoring binary."
        ub[vidx] = bval
    else
        # error in bound type
        @warn "Unknown bound type: $(btype). Input ignored."
    end

    return nothing
end