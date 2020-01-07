abstract type MPSSection end

struct MPSNoSection <: MPSSection end
struct MPSName <: MPSSection end
struct MPSObjsense <: MPSSection end
struct MPSObjName <: MPSSection end
struct MPSRows <: MPSSection end
struct MPSColumns <: MPSSection end
struct MPSRhs <: MPSSection end
struct MPSRanges <: MPSSection end
struct MPSBounds <: MPSSection end
struct MPSEndata <: MPSSection end

function MPSSection(sec)
    if sec == "NAME"
        return MPSName()
    elseif sec == "OBJSENSE"
        return MPSObjsense()
    elseif sec == "OBJNAME"
        return MPSObjName()
    elseif sec == "ROWS"
        return MPSRows()
    elseif sec == "COLUMNS"
        return MPSColumns()
    elseif sec == "RHS"
        return MPSRhs()
    elseif sec == "RANGES"
        return MPSRanges()
    elseif sec == "BOUNDS"
        return MPSBounds()
    elseif sec == "ENDATA"
        return MPSEndata()
    else
        return MPSNoSection()
    end
end


"""
    split_mps_line(s::String)

More efficient implementation than Julia's Base.split
"""
function split_mps_line(s::String)   
    buf = IOBuffer(s)
    S = String[]
    while !eof(buf)
        # Skip chars and read word
        skipchars(isspace, buf)
        n = position(buf)
        skipchars(!isspace, buf)
        m = position(buf)
        m > n || break
        push!(S, s[1+n:m])
    end
    
    return S
end

"""
    readmps!(m::Model{Tv}, filename)

Parse a free-MPS file into model `m`
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

    section = MPSNoSection()

    open(fname) do f
        nline = 0
        for ln in eachline(f)
            nline += 1

            # pass empty lines or comments
            if length(ln) == 0 || ln[1] == '*'
                continue
            end

            fields = split_mps_line(ln)

            # check for indicators
            if ln[1] != ' '
                # Section header
                section = MPSSection(fields[1])

                isa(section, MPSName) && (m.name = fields[2])
                continue
            end

            parseline!(section, m, fields, d, var2idx, con2idx, lb, ub)

            isa(section, MPSEndata) && break
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
    parseline!(::MPSSection, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Dummy function for now.
"""
function parseline!(::MPSSection, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    error()
end

"""
    parseline!(::MPSNoSection, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

This should never be called.
"""
function parseline!(::MPSNoSection, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    error()
end

"""
    parseline!(::MPSRows, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Parse a line from ROWS section in an MPS file.
"""
function parseline!(::MPSRows, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    length(fields) == 2 || error()

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
        # Error in the input
        error()
    end

    return nothing
end

"""
    parseline!(::MPSColumns, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Parse a line from COLUMNS section in an MPS file.
"""
function parseline!(::MPSColumns, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    # parse a ROWS line of the MPS file
    # `ln` is a string that contains the current line
    # `fields` contains the 6 fields of that line
    
    # current line should contain at least three fields
    length(fields) >= 3 || error("Line in COLUMNS section has less than three fields.")

    fields[2] == "'MARKER'" && (return nothing)

    # First field is empty, second field is variable's name
    cname = fields[1]
    if !haskey(var2idx, cname)
        # Create new variable
        vidx = add_variable!(m, String(cname), zero(Tv), zero(Tv), Tv(Inf))
        var2idx[cname] = vidx
        lb[vidx] = zero(Tv)
        ub[vidx] = Tv(Inf)
    end

    # Second and third fields are
    # the row's name, and the corresponding coefficient
    rname = fields[2]  # row name
    haskey(con2idx, rname) || error("Unknown row.")

    coeff = parse(Float64, fields[3])  # coefficient value
    if con2idx[rname].uuid == 0
        # objective coefficient
        set_obj_coeff!(m.pbdata_raw.vars[var2idx[cname]], coeff)
    else
        # constraint coefficient
        set_coeff!(m.pbdata_raw, var2idx[cname], con2idx[rname], coeff)
    end

    # optional other fields
    length(fields) >= 5 || (return nothing)

    rname = fields[4]  # row name
    haskey(con2idx, rname) || error("Unknown row.")
    coeff = parse(Float64, fields[5])  # coefficient value
    if con2idx[rname].uuid == 0
        # objective coefficient
        set_obj_coeff!(m.pbdata_raw.vars[var2idx[cname]], coeff)
    else
        # constraint coefficient
        set_coeff!(m.pbdata_raw, var2idx[cname], con2idx[rname], coeff)
    end

    return nothing
end

"""
    parseline!(::MPSRhs, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Parse a line from RHS section in an MPS file.
"""
function parseline!(::MPSRhs, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}

    # Line should contain at least 3 fields, and up to 5
    length(fields) >= 3 || error("Line is too short.")

    # First field is RHS name (may be empty)
    rhsname = fields[1]
    if d["rhs"] == ""
        # first time RHS is read
        d["rhs"] = rhsname
    elseif d["rhs"] != rhsname
        # other RHS, current line is ignored
        return nothing
    end

    # parse line
    rname = fields[2]
    haskey(con2idx, rname) || error()

    rval = parse(Float64, fields[3])
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

    if length(fields) < 5
        return nothing
    end

    # optional fields
    rname = fields[4]  # row name
    haskey(con2idx, rname) || error()
    rval = parse(Float64, fields[5])  # coefficient value
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

    return nothing
end

"""
    parseline!(::MPSRanges, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Parse a line from RANGES section in an MPS file.
"""
function parseline!(::MPSRanges, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    
    length(fields) >= 3 || error("Line too short")

    rngname = fields[1]
    if d["range"] == ""
        d["range"] = rngname
    elseif d["range"] != rngname
        # other range vector, ignore
        return nothing
    end

    # parse line
    rname = fields[2]
    rval = parse(Float64, fields[3])
    haskey(con2idx, rname) || error()
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

    if length(fields) >= 5
        rname = fields[4]
        rngval = parse(Float64, fields[5])
        haskey(con2idx, rname) || error()
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

    return nothing
end

"""
    parseline!(::MPSBounds, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

Parse a line from RANGES section in an MPS file.
"""
function parseline!(::MPSBounds, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    if length(fields) < 3 || (
            length(fields) < 4 && (
                fields[1] == "LO" || fields[1] == "UP" || fields[1] == "FX"
            )
        )
        # Missing fields
        error()
    end

    # first field should be bound type
    btype = fields[1]
    bname = fields[2]
    cname = fields[3]
    haskey(var2idx, cname) || error()

    # 
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
            lb[vidx] = 0.0
            ub[vidx] = 1.0
        else
            error("Unknown bound type: $(btype)")
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
        lb[vidx] = bval

    elseif btype == "UI"
        # 0 <= x <= b, x integer
        # Keep bounds but ignore integer requirement
        ub[vidx] = bval
    else
        # error in bound type
        error("Unknown bound type: $(btype)")
    end
    return nothing
end

"""
    parseline!(::MPSEndata, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub)

End of model. Do nothing
"""
function parseline!(::MPSEndata, m::Model{Tv}, fields, d, var2idx, con2idx, lb, ub) where{Tv<:Real}
    return nothing
end