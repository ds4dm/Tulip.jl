abstract type MPSSection end

struct MPSNoSection <: MPSSection end
struct MPSName <: MPSSection end
struct MPSObjsense <: MPSSection end
struct MPSRows <: MPSSection end
struct MPSColumns <: MPSSection end
struct MPSRhs <: MPSSection end
struct MPSRanges <: MPSSection end
struct MPSBounds <: MPSSection end
struct MPSEndata <: MPSSection end

function MPSSection(sec::String)
    if sec == "NAME"
        return MPSName()
    elseif sec == "OBJSENSE"
        return MPSObjsense()
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
    MPSData
"""
mutable struct MPSData{Tv<:Real}

    MPSFormat::Symbol  # :Free of :Fixed

    name::String  # Problem name
    objname::String  # Objective name
    rhsname::String  # name of RHS field
    boundsname::String  # name of BOUNDS field
    rangename::String  # nae of RANGES field

    sec::MPSSection

    varnames::Vector{String}  # idx -> name
    var2idx::Dict{String, Int}  # name -> index
    connames::Vector{String}  # idx -> name
    con2idx::Dict{String, Int}  # name -> index

    ncon::Int  # number of constraints (rows)
    nvar::Int  # Number of variables (columns)

    # Objective
    objsense::Symbol  # Objective sense
    c::Vector{Tv}  # Objective coefficients
    c0::Tv  # Objective offset

    # Coefficients in COO format
    aI::Vector{Int}
    aJ::Vector{Int}
    aV::Vector{Tv}

    # Bounds
    conbounds::Vector{Tuple{BoundType, Tv, Tv}}
    varbounds::Vector{Tuple{BoundType, Tv, Tv}}

    MPSData{Tv}() where{Tv<:Real}= new{Tv}(
        :Free,
        "", "", "", "", "", MPSNoSection(),
        String[], Dict{String, Int}(), String[], Dict{String, Int}(),
        0, 0,
        :Min, Tv[], zero(Tv),
        Int[], Int[], Tv[],
        Tuple{BoundType, Tv, Tv}[], Tuple{BoundType, Tv, Tv}[]
    )
end

"""
    split_mps_line(s::String)

More efficient implementation than Julia's Base.split
"""
function split_mps_line(s::String)   
    buf = IOBuffer(s)
    S = String[]
    sizehint!(S, 5)
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

"""
function readmps(fname::String)
    # MPS files have double float precision at most
    dat = MPSData{Float64}()

    mps = open(fname)
    try
        while !eof(mps)
            # Read line
            ln = readline(mps)

            if length(ln) == 0 || ln[1] == '*'
                # Line is empty or a comment
                continue
            end

            if !isspace(ln[1])
                # Section header
                fields = split_mps_line(ln)
                dat.sec = MPSSection(fields[1])
                isa(dat.sec, MPSEndata) && break
                if isa(dat.sec, MPSName)
                    dat.name = length(fields) >= 2 ? fields[2] : ""
                end
                continue
            end

            parseline!(dat.sec, dat, ln)
        end
    catch err
        close(mps)
        rethrow(err)
    end

    close(mps)

    # End of file reached
    isa(dat.sec, MPSEndata) || error("File ended without EOF flag.")

    return dat
end

function parseline!(::MPSObjsense, dat::MPSData, ln::String)
    # parse the name
    s = split_mps_line(ln)[1]
    if s == "MIN"
        dat.objsense = :Min
    elseif s == "MAX"
        dat.objsense = :Max
    else
        error("Unkown objective sense: $s")
    end
    return nothing
end

function parseline!(::MPSName, dat::MPSData, ln::String)
    # parse the name
    dat.objsense = split_mps_line(ln)[2]
    return nothing
end

"""
    parseline!(::MPSRows, dat, ln)

The current line is expected to have the form
```
 X ZZZZ 
```
where:
    * the first character is a space
    * `X` is either 'N' (objective row), `E` (equality constraint), `L` (less-or-equal), `G` (greater-or-equal)
    * `ZZZZ` cannot contain any space
"""
function parseline!(::MPSRows, dat::MPSData{Tv}, ln::String) where{Tv<:Real}
    buf = IOBuffer(ln)
    # Skip first spaces
    skipchars(isspace, buf)
    c = read(buf, Char)
    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    
    rname = ln[1+p1:p2]
    
    if c == 'N'
        if dat.objname == ""
            dat.objname = rname
        else
            # New objective row; ignore
        end
    elseif c == 'E' || c == 'L' || c == 'G'
        # constraint
        haskey(dat.con2idx, rname) && error("Dupplicate row name $rname")
        dat.ncon += 1

        # Record constraint name
        push!(dat.connames, rname)
        dat.con2idx[rname] = dat.ncon

        # Create default bounds
        if c == 'E'
            push!(dat.conbounds, (TLP_FX, zero(Tv), zero(Tv)))
        elseif c == 'L'
            push!(dat.conbounds, (TLP_UP, Tv(-Inf), zero(Tv)))
        elseif c == 'G'
            push!(dat.conbounds, (TLP_LO, zero(Tv), Tv(Inf)))
        end
    else
        # Error in the input
        error("Unknown row type: $c")
    end

    return nothing
end

"""
    parseline!(::MPSColumns, dat, ln)

The current line is expected to have the form
```
```
"""
function parseline!(::MPSColumns, dat::MPSData{Tv}, ln::String) where{Tv<:Real}
    fields = split_mps_line(ln)

    length(fields) >= 3 || error("Incomplete line")
    vname = fields[1]
    rname = fields[2]

    # Ignore Intger markers
    rname == "'MARKER'" && return nothing

    coeff = parse(Tv, fields[3])

    if !haskey(dat.var2idx, vname)
        # Create new variable
        dat.nvar += 1
        dat.var2idx[vname] = dat.nvar
        push!(dat.varnames, vname)

        # Create default bounds and zero objective coeff
        push!(dat.varbounds, (TLP_LO, zero(Tv), Tv(Inf)))

        push!(dat.c, zero(Tv))

        j = dat.nvar
    else
        # Variable already exists; record its index
        j = dat.var2idx[vname]
    end

    # Add coefficient
    if dat.objname == rname
        dat.c[j] = coeff
    elseif haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]
        push!(dat.aI, i)
        push!(dat.aJ, j)
        push!(dat.aV, coeff)
    else
        # Ignore input
    end

    length(fields) >= 5 || return nothing  # end of line reached

    # Read fields 4 and 5
    rname = fields[4]
    coeff = parse(Tv, fields[5])

    # Add coefficient
    if dat.objname == rname
        dat.c[j] = coeff
    elseif haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]
        push!(dat.aI, i)
        push!(dat.aJ, j)
        push!(dat.aV, coeff)
    else
        # Unknown row
        # error("Unkown row $rname")
    end

    return nothing
end

"""
    parseline!(::MPSRhs, dat, ln)

The current line is expected to have the form
```
```
"""
function parseline!(::MPSRhs, dat::MPSData{Tv}, ln::String) where{Tv<:Real}
    buf = IOBuffer(ln)
    # Skip first spaces
    skipchars(isspace, buf)

    # Read RHS name field
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")

    rhsname = ln[1+p1:p2]

    if dat.rhsname == ""
        dat.rhsname = rhsname
    elseif dat.rhsname != rhsname
        # Other RHS, current line is ignored
        return nothing
    end

    # Read fields 2 and 3
    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")
    rname = ln[1+p1:p2]

    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    rval = parse(Tv, ln[1+p1:p2])

    if dat.objname == rname
        # Objective offset
        i = 0
        dat.c0 = -rval
    elseif haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]

        # update constraint bounds
        bt, lb, ub = dat.conbounds[i]
        if bt == TLP_UP
            # a'x <= b
            dat.conbounds[i] = (bt, Tv(-Inf), rval)
        elseif bt == TLP_LO
            # a'x >= b
            dat.conbounds[i] = (bt, rval, Tv(Inf))
        elseif bt == TLP_FX
            # a'x == b
            dat.conbounds[i] = (bt, rval, rval)
        else
            error("Got RHS term for row $rname of type $bt")
        end
    else
        error("Unkown row $rname")
    end


    skipchars(isspace, buf)
    eof(buf) && return nothing  # end of line reached

    # Read fields 4 and 5
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")
    rname = ln[1+p1:p2]

    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    rval = parse(Tv, ln[1+p1:p2])

    if dat.objname == rname
        # Objective offset
        i = 0
        dat.c0 = -rval
    elseif haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]

        # update constraint bounds
        bt, lb, ub = dat.conbounds[i]
        if bt == TLP_UP
            # a'x <= b
            dat.conbounds[i] = (bt, Tv(-Inf), rval)
        elseif bt == TLP_LO
            # a'x >= b
            dat.conbounds[i] = (bt, rval, Tv(Inf))
        elseif bt == TLP_FX
            # a'x == b
            dat.conbounds[i] = (bt, rval, rval)
        else
            error("Got RHS term for row $rname of type $bt")
        end
    else
        error("Unkown row $rname")
    end

    return nothing
end

"""
    parseline!(::MPSRanges, dat, ln)

The current line is expected to have the form
```
```
"""
function parseline!(::MPSRanges, dat::MPSData{Tv}, ln::String) where{Tv<:Real}
    buf = IOBuffer(ln)
    # Skip first spaces
    skipchars(isspace, buf)

    # Read RANGE name field
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")

    rgname = ln[1+p1:p2]

    if dat.rangename == ""
        dat.rangename = rgname
    elseif dat.rangename != rgname
        # Other RHS, current line is ignored
        return nothing
    end

    # Read fields 2 and 3
    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")
    rname = ln[1+p1:p2]

    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    rval = parse(Tv, ln[1+p1:p2])

    if haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]

        # update constraint bounds
        bt, lb, ub = dat.conbounds[i]
        if bt == TLP_UP
            # `-Inf <= a'x <= u` becomes `u - |r| <= a'x <= u`
            dat.conbounds[i] = (TLP_RG, ub - abs(rval), ub)
        elseif bt == TLP_LO
            # `l <= a'x <= Inf` becomes `l <= a'x <= l + |r|`
            dat.conbounds[i] = (TLP_RG, lb, lb + abs(rval))
        elseif bt == TLP_FX && rval >= zero(Tv)
            # `a'x = b` becomes `b <= a'x <= b + |r|`
            dat.conbounds[i] = (TLP_RG, lb, lb + rval)
        elseif bt == TLP_FX && rval < zero(Tv)
            # `a'x = b` becomes `b - |r| <= a'x <= b`
            dat.conbounds[i] = (TLP_RG, lb + rval, lb)
        else
            error("Got RANGE val $rval for row $rname of type $bt")
        end
    else
        error("Unkown row $rname")
    end

    skipchars(isspace, buf)
    eof(buf) && return nothing  # end of line reached

    # Read fields 4 and 5
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    eof(buf) && error("Incomplete line")
    rname = ln[1+p1:p2]

    skipchars(isspace, buf)
    p1 = position(buf)
    skipchars(!isspace, buf)
    p2 = position(buf)
    rval = parse(Tv, ln[1+p1:p2])

    if haskey(dat.con2idx, rname)
        i = dat.con2idx[rname]

        # update constraint bounds
        bt, lb, ub = dat.conbounds[i]
        if bt == TLP_UP
            # `-Inf <= a'x <= u` becomes `u - |r| <= a'x <= u`
            dat.conbounds[i] = (TLP_RG, ub - abs(rval), ub)
        elseif bt == TLP_LO
            # `l <= a'x <= Inf` becomes `l <= a'x <= l + |r|`
            dat.conbounds[i] = (TLP_RG, lb, lb + abs(rval))
        elseif bt == TLP_FX && rval >= zero(Tv)
            # `a'x = b` becomes `b <= a'x <= b + |r|`
            dat.conbounds[i] = (TLP_RG, lb, lb + rval)
        elseif bt == TLP_FX && rval < zero(Tv)
            # `a'x = b` becomes `b - |r| <= a'x <= b`
            dat.conbounds[i] = (TLP_RG, lb + rval, lb)
        else
            error("Got RANGE val $rval for row $rname of type $bt")
        end
    else
        error("Unkown row $rname")
    end

    return nothing
end

"""
    parseline!(::MPSBounds, dat, ln)

The current line is expected to have the form
```
```
"""
function parseline!(::MPSBounds, dat::MPSData{Tv}, ln::String) where{Tv<:Real}
    fields = split_mps_line(ln)

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
    bdname = fields[2]
    colname = fields[3]

    if dat.boundsname == ""
        dat.boundsname = bdname
    elseif dat.boundsname != bdname
        # Other BOUNDS, current line is ignored
        return nothing
    end

    haskey(dat.var2idx, colname) || error()

    # Get column index
    j = dat.var2idx[colname]
    bt, lb, ub = dat.varbounds[j]

    if length(fields) == 3
        if btype == "FR"
            # -∞ < x < ∞
            dat.varbounds[j] = (TLP_FR, Tv(-Inf), Tv(Inf))
        elseif btype == "MI"
            # -∞ < x <= ub
            if isfinite(ub)
                dat.varbounds[j] = (TLP_UP, Tv(-Inf), ub)
            else
                dat.varbounds[j] = (TLP_FR, Tv(-Inf), Tv(Inf))
            end
        elseif btype == "PL"
            # lb <= x < ∞
            if isfinite(lb)
                dat.varbounds[j] = (TLP_LO, lb, Tv(Inf))
            else
                dat.varbounds[j] = (TLP_FR, Tv(-Inf), Tv(Inf))
            end
        elseif btype == "BV"
            # x ∈ {0, 1}
            # Record bounds 0 <= x <= 1 but ignore binary requirement
            dat.varbounds[j] = (TLP_RG, zero(Tv), one(Tv))
        else
            error("Unknown bound type: $(btype)")
        end

        return nothing
    end

    bval = parse(Tv, fields[4])
    
    if btype == "LO" || btype == "LI"
        # b <= x <= ub
        # Ignore integer requirement
        if isfinite(ub)
            dat.varbounds[j] = (TLP_RG, bval, ub)
        else
            dat.varbounds[j] = (TLP_LO, bval, Tv(Inf))
        end
    
    elseif btype == "UP" || btype == "UI"
        # lb <= x <= b
        # Ignore integer requirement
        if isfinite(lb)
            dat.varbounds[j] = (TLP_RG, lb, bval)
        else
            dat.varbounds[j] = (TLP_UP, -Tv(Inf), bval)
        end
    
    elseif btype == "FX"
        # x == b
        dat.varbounds[j] = (TLP_FX, bval, bval)
    
    else
        # error in bound type
        error("Unknown bound type: $(btype)")
    end
    return nothing
end