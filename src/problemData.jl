using SparseArrays

mutable struct RowOrCol{Tv}
    nzind::Vector{Int}
    nzval::Vector{Tv}
end

const Row = RowOrCol
const Col = RowOrCol

"""
    ProblemData{Tv}

Data structure for storing problem data in precision `Tv`.

The LP is represented in canonical form

```math
\\begin{array}{rl}
    \\displaystyle \\min_{x} \\ \\ \\ & c^{T} x + c_{0} \\\\
    s.t. \\ \\ \\ & l_{r} \\leq A x \\leq u_{r} \\\\
    & l_{c} \\leq x \\leq u_{c}
\\end{array}
```
"""
mutable struct ProblemData{Tv}

    name::String

    # Dimensions
    ncon::Int  # Number of rows
    nvar::Int  # Number of columns (i.e. variables)

    # Objective
    # TODO: objective sense
    objsense::Bool  # true is min, false is max
    obj::Vector{Tv}
    obj0::Tv  # Constant objective offset

    # Constraint matrix
    # We store both rows and columns. It is redundant but simplifies access.
    # TODO: put this in its own data structure? (would allow more flexibility in modelling)
    arows::Vector{Row{Tv}}
    acols::Vector{Col{Tv}}

    # TODO: Data structures for QP
    # qrows
    # qcols

    # Bounds
    lcon::Vector{Tv}
    ucon::Vector{Tv}
    lvar::Vector{Tv}
    uvar::Vector{Tv}

    # Names
    con_names::Vector{String}
    var_names::Vector{String}

    # Only allow empty problems to be instantiated for now
    ProblemData{Tv}(pbname::String="") where {Tv} = new{Tv}(
        pbname, 0, 0,
        true, Tv[], zero(Tv),
        Row{Tv}[], Col{Tv}[],
        Tv[], Tv[], Tv[], Tv[],
        String[], String[]
    )
end

import Base.empty!

function Base.empty!(pb::ProblemData{Tv}) where{Tv}
    
    pb.name = ""
    
    pb.ncon = 0
    pb.nvar = 0
    
    pb.objsense = true
    pb.obj = Tv[]
    pb.obj0 = zero(Tv)
    
    pb.arows = Row{Tv}[]
    pb.acols = Col{Tv}[]
    
    pb.lcon = Tv[]
    pb.ucon = Tv[]
    pb.lvar = Tv[]
    pb.uvar = Tv[]

    pb.con_names = String[]
    pb.var_names = String[]
    
    return pb
end

#=
    TODO:
    * Creation
        *[x] Add single constraint
        *[ ] Add multiple constraints
        *[x] Add single variable
        *[ ] Add multiple variables
        *[x] Load entire problem
    * Modification
        *[x] Empty model
        *[x] Delete single constraint
        *[x] Delete multiple constraints (fallback)
        *[x] Delete single variable
        *[x] Delete multiple variables (fallback)
        *[x] Change single coefficient
        *[ ] Change multiple coefficients
    * Attributes
        *[ ] Query model attributes
            *[ ] MOI-supported attributes
            *[ ] Other attributes
=#

# =============================
#     Problem creation
# =============================

"""
    add_constraint!(pb, rind, rval, l, u; [name, issorted])

Add one linear constraint to the problem.

# Arguments
* `pb::ProblemData{Tv}`: the problem to which the new row is added
* `rind::Vector{Int}`: column indices in the new row
* `rval::Vector{Tv}`: non-zero values in the new row
* `l::Tv`
* `u::Tv`
* `name::String`: row name (defaults to `""`)
* `issorted::Bool`: indicates whether the row indices are already issorted.
"""
function add_constraint!(pb::ProblemData{Tv},
    rind::Vector{Int}, rval::Vector{Tv},
    l::Tv, u::Tv,
    name::String="";
    issorted::Bool=false
)::Int where{Tv}
    # Sanity checks
    nz = length(rind)
    nz == length(rval) || throw(DimensionMismatch(
        "Cannot add a row with $nz indices but $(length(rval)) non-zeros"
    ))

    # Increment row counter
    pb.ncon += 1
    push!(pb.con_names, name)
    push!(pb.lcon, l)
    push!(pb.ucon, u)

    if nz == 0
        # emtpy row
        push!(pb.arows, Row{Tv}(Int[], Tv[]))
        return pb.ncon
    end

    # TODO: check coefficients values, e.g.
    #   * all coefficients are finite
    #   * remove hard zeros
    #   * combine dupplicate indices
    #   * check if indices are already issorted
    
    # Create new row
    if issorted
        row = Row{Tv}(copy(rind), copy(rval))
    else
        # Sort indices first
        p = sortperm(rind)
        row = Row{Tv}(rind[p], rval[p])
    end
    push!(pb.arows, row)

    # Update column coefficients
    for (j, v) in zip(rind, rval)
        if !iszero(v)
            push!(pb.acols[j].nzind, pb.ncon)
            push!(pb.acols[j].nzval, v)
        end
    end

    # Done
    return pb.ncon
end

"""
    add_variable!(pb, cind, cval, obj, l, u, [name])

Add one variable to the problem.

# Arguments
* `pb::ProblemData{Tv}`: the problem to which the new column is added
* `cind::Vector{Int}`: row indices in the new column
* `cval::Vector{Tv}`: non-zero values in the new column
* `obj::Tv`: objective coefficient
* `l::Tv`: column lower bound
* `u::Tv`: column upper bound
* `name::String`: column name (defaults to `""`)
* `issorted::Bool`: indicates whether the column indices are already issorted.
"""
function add_variable!(pb::ProblemData{Tv},
    cind::Vector{Int}, cval::Vector{Tv},
    obj::Tv, l::Tv, u::Tv,
    name::String="";
    issorted::Bool=false
)::Int where{Tv}
    # Sanity checks
    nz = length(cind)
    nz == length(cval) || throw(DimensionMismatch(
        "Cannot add a column with $nz indices but $(length(cval)) non-zeros"
    ))

    # Increment column counter
    pb.nvar += 1
    push!(pb.var_names, name)
    push!(pb.lvar, l)
    push!(pb.uvar, u)
    push!(pb.obj, obj)

    if nz == 0
        push!(pb.acols, Col{Tv}(Int[], Tv[]))
        return pb.nvar  # empty column
    end

    # TODO: check coefficients

    # Create a new column
    if issorted
        col = Col{Tv}(copy(cind), copy(cval))
    else
        # Sort indices
        p = sortperm(cind)
        col = Col{Tv}(cind[p], cind[p])
    end
    push!(pb.acols, col)

    # Update row coefficients
    for (i, v) in zip(cind, cval)
        if !iszero(v)
            push!(pb.arows[i].nzind, pb.nvar)
            push!(pb.arows[i].nzval, v)
        end
    end

    # Done
    return pb.nvar
end

"""
    load_problem!(pb, )

Load entire problem.
"""
function load_problem!(pb::ProblemData{Tv},
    name::String,
    objsense::Bool, obj::Vector{Tv}, obj0::Tv,
    A::SparseMatrixCSC,
    lcon::Vector{Tv}, ucon::Vector{Tv},
    lvar::Vector{Tv}, uvar::Vector{Tv},
    con_names::Vector{String}, var_names::Vector{String}
) where{Tv}
    empty!(pb)

    # Sanity checks
    ncon, nvar = size(A)
    ncon == length(lcon) || error("")
    ncon == length(ucon) || error("")
    ncon == length(con_names) || error("")
    nvar == length(obj)
    isfinite(obj0) || error("Objective offset $obj0 is not finite")
    nvar == length(lvar) || error("")
    nvar == length(uvar) || error("")

    # Copy data
    pb.name = name
    pb.ncon = ncon
    pb.nvar = nvar
    pb.objsense = objsense
    pb.obj = copy(obj)
    pb.obj0 = obj0
    pb.lcon = copy(lcon)
    pb.ucon = copy(ucon)
    pb.lvar = copy(lvar)
    pb.uvar = copy(uvar)
    pb.con_names = copy(con_names)
    pb.var_names = copy(var_names)

    # Load coefficients
    pb.acols = Vector{Col{Tv}}(undef, nvar)
    pb.arows = Vector{Row{Tv}}(undef, ncon)
    for j in 1:nvar
        col = A[:, j]
        pb.acols[j] = Col{Tv}(col.nzind, col.nzval)
    end

    At = sparse(A')
    for i in 1:ncon
        row = At[:, i]
        pb.arows[i] = Row{Tv}(row.nzind, row.nzval)
    end

    return pb
end

# =============================
#     Problem modification
# =============================

"""
    delete_constraint!(pb::ProblemData, rind::Int)

Delete a single constraint from problem `pb`.
"""
function delete_constraint!(pb::ProblemData{Tv}, rind::Int) where{Tv}
    # Sanity checks
    1 <= rind <= pb.ncon || error("Invalid row index $rind")

    # Delete row name and bounds
    deleteat!(pb.con_names, rind)
    deleteat!(pb.lcon, rind)
    deleteat!(pb.ucon, rind)
    
    # Update columns
    for (j, col) in enumerate(pb.acols)
        # Search for row in that column
        rg = searchsorted(col.nzind, rind)
        if rg.start > length(col.nzind)
            # Nothing to do
            continue
        else
            if col.nzind[rg.start] == rind
                # Delete row from column
                deleteat!(col.nzind, rg.start)
                deleteat!(col.nzval, rg.start)
            end

            # Decrement subsequent row indices
            col.nzind[rg.start:end] .-= 1
        end
    end 

    # Delete row
    deleteat!(pb.arows, rind)

    # Update row counter
    pb.ncon -= 1
    return nothing
end

"""
    delete_constraints!(pb::ProblemData, rinds)

Delete rows in collection `rind` from problem `pb`.

# Arguments
* `pb::ProblemData`
* `rinds`: collection of row indices to be removed
"""
function delete_constraints!(pb::ProblemData{Tv}, rinds) where{Tv}
    # TODO: don't use fallback 
    for i in rinds
        delete_constraint!(pb, i)
    end
    return nothing
end

"""
    delete_variable!(pb, cind)

Delete a single column from problem `pb`.
"""
function delete_variable!(pb::ProblemData{Tv}, cind::Int) where{Tv}
    # Sanity checks
    1 <= cind <= pb.nvar || error("Invalid column index $cind")

    # Delete column name, objective and bounds
    deleteat!(pb.var_names, cind)
    deleteat!(pb.obj,  cind)
    deleteat!(pb.lvar, cind)
    deleteat!(pb.uvar, cind)
    
    # Update rows
    for (i, row) in enumerate(pb.arows)
        # Search for column in that row
        rg = searchsorted(row.nzind, cind)
        if rg.start > length(row.nzind)
            # Nothing to do
            continue
        else
            if row.nzind[rg.start] == cind
                # Column appears in row
                deleteat!(row.nzind, rg.start)
                deleteat!(row.nzval, rg.start)
            end

            # Decrement subsequent column indices
            row.nzind[rg.start:end] .-= 1
        end
    end

    # Delete column
    deleteat!(pb.acols, cind)

    # Update column counter
    pb.nvar -= 1
    return nothing
end

"""
    delete_variables!(pb::ProblemData, cinds)

Delete a collection of columns from problem `pb`.

# Arguments
* `pb::ProblemData`
* `cinds`: collection of row indices to be removed
"""
function delete_variables!(pb::ProblemData{Tv}, cinds) where{Tv}
    # TODO: don't use fallback 
    for j in cinds
        delete_variable!(pb, j)
    end
    return nothing
end

"""
    set_coefficient!(pb, i, j, v)

Set the coefficient `(i, j)` to value `v`.

# Arguments
* `pb::ProblemData{Tv}`: the problem whose coefficient
* `i::Int`: row index
* `j::Int`: column index
* `v::Tv`: coefficient value
"""
function set_coefficient!(pb::ProblemData{Tv}, i::Int, j::Int, v::Tv) where{Tv}
    # Sanity checks
    1 <= i <= pb.ncon && 1 <= j <= pb.nvar || error(
        "Cannot access coeff $((i, j)) in a model of size ($(pb.ncon), $(pb.nvar))"
    )

    # Update row and column
    _set_coefficient!(pb.arows[i], j, v)
    _set_coefficient!(pb.acols[j], i, v)

    return nothing
end

"""
    _set_coefficient!(roc::RowOrCol{Tv}, ind::Int, v::Tv)

Set coefficient to value `v`.
"""
function _set_coefficient!(roc::RowOrCol{Tv}, ind::Int, v::Tv) where{Tv}
    # Check if index already exists
    k = searchsortedfirst(roc.nzind, ind)

    if (1 <= k <= length(roc.nzind)) && roc.nzind[k] == ind
        # This coefficient was a non-zero before
        if iszero(v)
            deleteat!(roc.nzind, k)
            deleteat!(roc.nzval, k)
        else
            roc.nzval[k] = v
        end
    else
        # Only add coeff if non-zero
        if !iszero(v)
            insert!(roc.nzind, k, ind)
            insert!(roc.nzval, k, v)
        end
    end

    return nothing
end

# =============================
#     Problem queries
# =============================

# TODO