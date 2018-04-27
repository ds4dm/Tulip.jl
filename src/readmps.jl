"""
    readmps(filename)

    Parse an MPS file, and return a linear program in standard form.
"""
function readmps(file_name)
    # parse an MPS file
    section = ""
    d = Dict()

    row2idx = Dict()  # name <-> index correspondence for rows
    col2idx = Dict()  # name <-> index correspondence for columns

    senses = []  # sense of rows
     
    coeffs_row = Array{Int, 1}(0)
    coeffs_col = Array{Int, 1}(0)
    coeffs_val = Array{Float64, 1}(0)

    obj_col = Array{Int, 1}(0)
    obj_val = Array{Float64, 1}(0)

    rhs_row = Array{Int, 1}(0)
    rhs_val = Array{Float64, 1}(0)

    lb_col = Array{Int, 1}(0)
    lb_val = Array{Float64, 1}(0)
    ub_col = Array{Int, 1}(0)
    ub_val = Array{Float64, 1}(0)

    ranges_row = Array{Int, 1}(0)
    ranges_val = Array{Float64, 1}(0)

    d["nrows"] = 0
    d["ncols"] = 0
    d["rhs"] = ""
    d["bound"] = ""
    d["range"] = ""

    open(file_name) do f
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
                        d["name"] = fields[2]
                    else
                        #  empty name
                        d["name"] = ""
                    end
                    continue
                end
                # update section value and pass to next line
                section = fields[1]
                continue
            end

            if section == "ROWS"
                parsemps_rows!(ln, fields, d, row2idx, senses)

            elseif section == "COLUMNS"
                parsemps_columns!(
                    ln, fields, d, col2idx, row2idx,
                    coeffs_row, coeffs_col, coeffs_val, obj_col, obj_val
                )

            elseif section == "RHS"
                parsemps_rhs!(ln, fields, d, row2idx, rhs_row, rhs_val)

            elseif section == "RANGES"
                parsemps_ranges!(ln, fields, d, row2idx, ranges_row, ranges_val)

            elseif section == "BOUNDS"
                parsemps_bounds!(ln, fields, d, col2idx, lb_col, lb_val, ub_col, ub_val)

            elseif section == "ENDATA"
                # end of file
            else
                warn("UNEXPECTED FORMAT PARSING\t$(ln)")
            end

        end

        # End of file reached
    end

    # extract relevant info
    ranges = sparsevec(ranges_row, ranges_val)

    # transform to standard form
    m, n, obj, rhs, coeffs, lb, ub = convert_to_standard_form(
        d, 
        obj_col, obj_val,
        rhs_row, rhs_val, senses,
        coeffs_row, coeffs_col, coeffs_val,
        lb_col, lb_val, ub_col, ub_val,
        ranges_row, ranges_val
    )
    
    return m, n, obj, rhs, coeffs, lb, ub, ranges
end


function convert_to_standard_form(
    d, 
    obj_col, obj_val,  # objective
    rhs_row, rhs_val, senses,  # right-hand side and senses
    coeffs_row, coeffs_col, coeffs_val,  # coefficients
    lb_col, lb_val, ub_col, ub_val,  # lower and uper bounds
    ranges_row, ranges_val
    )
    
    m = d["nrows"]  # number of constraints in original formulation
    n = d["ncols"]  # number of variables in original formulation

    # I. Tranform constraints to equality constraints
    for j=1:m
        if senses[j] == "L"
            # add slack variable
            n += 1
            push!(coeffs_row, j)
            push!(coeffs_col, n)
            push!(coeffs_val, 1.0)
            senses[j] = "E"

        elseif senses[j] == "G"
            # add surplus variable
            n += 1
            push!(coeffs_row, j)
            push!(coeffs_col, n)
            push!(coeffs_val, -1.0)
            senses[j] = "E"
        end
    end

    coeffs = sparse(coeffs_row, coeffs_col, coeffs_val, m, n)
    obj = Array(sparsevec(obj_col, obj_val, n))
    rhs = Array(sparsevec(rhs_row, rhs_val, m))


    # II. Form bounds
    lb = zeros(n)  # lower bound
    ub = Inf * ones(n)  # upper bound
    for i in lb_col
        lb[i] = lb_val[i] 
    end
    for i in ub_col
        ub[i] = ub_val[i]
    end

    return m, n, obj, rhs, coeffs, lb, ub
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
        
        if length(f_) == 0 && length(s) == 0
            # f_ is empty field at the right-end
            # do nothing
            continue
        # elseif contains(f_, "\$")
        #     # f contains the $ sign: remove everything after $
        #     s = [split(f_, '$')[1]]
        #     # rest of line is treated as comment
        #     continue
        else
            # f_ is empty but there is a non-empty field later
            unshift!(s, f_)
            continue
        end
    end

    return s
end


function parsemps_rows!(ln, fields, d, row2idx, senses)
    # parse a ROWS line of the MPS file
    # `ln` is a string that contains the current line
    # `fields` is an array such that fields == split(ln)
    
    # current line should contain exactly 2 fields
    if length(fields) != 2
        # input error, current line is ignored
        warn("INPUT ERROR:\t$(ln)")
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
        row2idx[fields[2]] = 0

    elseif fields[1] == "E" || fields[1] == "L" || fields[1] == "G"
        # constraint
        d["nrows"] += 1
        row2idx[fields[2]] = d["nrows"]  # name
        push!(senses, fields[1])  # sense (i.e., E, L, or G)

    else
        # input error, current line is ignored
        warn("INPUT ERROR:\t$(ln)")
        return nothing
    end

    return nothing
end


function parsemps_columns!(ln, fields, d, col2idx, row2idx,
    coeffs_row, coeffs_col, coeffs_val, obj_col, obj_val
)
    # parse a ROWS line of the MPS file
    # `ln` is a string that contains the current line
    # `fields` contains the 6 fields of that line
    
    # current line should contain at least three fields
    if length(fields) < 4
        warn("COLUMN LINE TOO SHORT|$(ln)")
        return nothing
    end

    # First field is empty, second field is variable's name
    cname = fields[2]
    if !haskey(col2idx, cname)
        d["ncols"] += 1
        col2idx[cname] = d["ncols"]  # Add column name to the pool
    end

    # Second and third fields are
    # the row's name, and the corresponding coefficient
    rname = fields[3]  # row name
    if !haskey(row2idx, rname)
        # current row not in list of rows, current line ignored
        warn("COL ERROR UNKNOWN ROW $(rname)|$(ln)")
        return nothing
    end
    coeff = parse(Float64, fields[4])  # coefficient value
    if row2idx[rname] == 0
        # objective coefficient
        push!(obj_col, col2idx[cname])
        push!(obj_val, coeff)
    else
        # constraint coefficient
        push!(coeffs_col, col2idx[cname])
        push!(coeffs_row, row2idx[rname])
        push!(coeffs_val, coeff)
    end

    # optional other fields
    if length(fields) >= 6
        rname = fields[5]  # row name
        if !haskey(row2idx, rname)
            # current row not in list of rows, current line ignored
            warn("COL ERROR UNKNOWN ROW $(rname)|$(ln)")
            return nothing
        end
        coeff = parse(Float64, fields[6])  # coefficient value
        if row2idx[rname] == 0
            # objective coefficient
            push!(obj_col, col2idx[cname])
            push!(obj_val, coeff)
        else
            # constraint coefficient
            push!(coeffs_col, col2idx[cname])
            push!(coeffs_row, row2idx[rname])
            push!(coeffs_val, coeff)
        end
    end

    return nothing
end


function parsemps_rhs!(ln, fields, d, row2idx, rhs_row, rhs_val)
    # parse a line of the RHS section
    # `ln` is the current line
    # `fields` contains the fields of that line


    # Check for too short lines
    if length(fields) < 4
        # input error, current line is ignored
        warn("RHS LINE TOO SHORT|$(ln)")
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
    if !haskey(row2idx, rname)
        # current row not in list of rows, current line ignored
        warn("RHS ERROR UNKNOWN ROW $(rname)|$(ln)")
        return nothing
    end
    rval = parse(Float64, fields[4])
    # update index and value
    if row2idx[rname] == 0
        d["obj_offset"] = -rval
    else
        push!(rhs_row, row2idx[rname])
        push!(rhs_val, rval)
    end

    # optional fields
    if length(fields) >= 6
        rname = fields[5]  # row name
        if !haskey(row2idx, rname)
            # current row not in list of rows, current line ignored
            warn("RHS ERROR UNKNOWN ROW $(rname)|$(ln)")
            return nothing
        end
        rval = parse(Float64, fields[6])  # coefficient value
        # update index and value
        if row2idx[rname] == 0
            d["obj_offset"] = -rval
        else
            push!(rhs_row, row2idx[rname])
            push!(rhs_val, rval)
        end
    end

    return nothing
end


function parsemps_ranges!(ln, fields, d, row2idx, ranges_row, ranges_val)
    if length(fields) < 4
        # input error, current line is ignored
        warn("RNG LINE TOO SHORT|$(ln)")
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
    rngval = parse(Float64, fields[4])
    if !haskey(row2idx, rname)
        # unknown row
        warn("RNG ERROR UNKNOWN ROW $(rname)|$(ln)")
        return nothing
    end
    push!(ranges_row, row2idx[rname])
    push!(ranges_val, rngval)

    if length(fields) >=6
        rname = fields[5]
        rngval = parse(Float64, fields[6])
        if !haskey(row2idx, rname)
            # unknown row
            warn("RNG ERROR UNKNOWN ROW $(rname)|$(ln)")
            return nothing
        end
        push!(ranges_row, row2idx[rname])
        push!(ranges_val, rngval)
    end

    return nothing
end


function parsemps_bounds!(ln, fields, d, col2idx, lb_col, lb_val, ub_col, ub_val)
    # parse a line of the BOUNDS section

    if (length(fields) < 3
        || ((length(fields) < 4) 
            && (fields[1] == "LO" || fields[1] == "UP" || fields[1] == "FX"))
    )
        warn("BND ERROR LINE TOO SHORT|$(ln)")
        return nothing
    end

    # first field should be bound type
    btype = fields[1]
    bname = fields[2]
    cname = fields[3]
    if !haskey(col2idx, cname)
        warn("BND ERROR UNKNOWN COL $(rname)|$(ln)")
        return nothing
    end

    if btype == "LO"
        push!(lb_col, col2idx[cname])
        push!(lb_val, parse(Float64, fields[4]))
    elseif btype == "UP"
        push!(ub_col, col2idx[cname])
        push!(ub_val, parse(Float64, fields[4]))
    elseif btype == "FR"
        push!(lb_col, col2idx[cname])
        push!(lb_val, -Inf)
        push!(ub_col, col2idx[cname])
        push!(ub_val, Inf)
    elseif btype == "FX"
        push!(lb_col, col2idx[cname])
        push!(lb_val, parse(Float64, fields[4]))
        push!(ub_col, col2idx[cname])
        push!(ub_val, parse(Float64, fields[4]))
    elseif btype == "MI"
        push!(lb_col, col2idx[cname])
        push!(lb_val, -Inf)
    elseif btype == "PL"
        push!(ub_col, col2idx[cname])
        push!(ub_val, Inf)
    else
        # error in bound type
        warn("BND ERROR WRONG BOUND TYPE $(btype)")
    end

    return nothing
end