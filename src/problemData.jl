using SparseArrays

mutable struct RowOrCol{T}
    nzind::Vector{Int}
    nzval::Vector{T}
end

const Row = RowOrCol
const Col = RowOrCol

"""
    ProblemData{T}

Data structure for storing problem data in precision `T`.

The LP is represented in canonical form

```math
\\begin{array}{rl}
    \\displaystyle \\min_{x} \\ \\ \\ & c^{T} x + c_{0} \\\\
    s.t. \\ \\ \\ & l_{r} \\leq A x \\leq u_{r} \\\\
    & l_{c} \\leq x \\leq u_{c}
\\end{array}
```
"""
mutable struct ProblemData{T}

    name::String

    # Dimensions
    ncon::Int  # Number of rows
    nvar::Int  # Number of columns (i.e. variables)

    # Objective
    # TODO: objective sense
    objsense::Bool  # true is min, false is max
    obj::Vector{T}
    obj0::T  # Constant objective offset

    # Constraint matrix
    # We store both rows and columns. It is redundant but simplifies access.
    # TODO: put this in its own data structure? (would allow more flexibility in modelling)
    arows::Vector{Row{T}}
    acols::Vector{Col{T}}

    # TODO: Data structures for QP
    # qrows
    # qcols

    # Bounds
    lcon::Vector{T}
    ucon::Vector{T}
    lvar::Vector{T}
    uvar::Vector{T}

    # Names
    con_names::Vector{String}
    var_names::Vector{String}

    # Only allow empty problems to be instantiated for now
    ProblemData{T}(pbname::String="") where {T} = new{T}(
        pbname, 0, 0,
        true, T[], zero(T),
        Row{T}[], Col{T}[],
        T[], T[], T[], T[],
        String[], String[]
    )
end

import Base.empty!

function Base.empty!(pb::ProblemData{T}) where{T}

    pb.name = ""

    pb.ncon = 0
    pb.nvar = 0

    pb.objsense = true
    pb.obj = T[]
    pb.obj0 = zero(T)

    pb.arows = Row{T}[]
    pb.acols = Col{T}[]

    pb.lcon = T[]
    pb.ucon = T[]
    pb.lvar = T[]
    pb.uvar = T[]

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
* `pb::ProblemData{T}`: the problem to which the new row is added
* `rind::Vector{Int}`: column indices in the new row
* `rval::Vector{T}`: non-zero values in the new row
* `l::T`
* `u::T`
* `name::String`: row name (defaults to `""`)
* `issorted::Bool`: indicates whether the row indices are already issorted.
"""
function add_constraint!(pb::ProblemData{T},
    rind::Vector{Int}, rval::Vector{T},
    l::T, u::T,
    name::String="";
    issorted::Bool=false
)::Int where{T}
    # Sanity checks
    nz = length(rind)
    nz == length(rval) || throw(DimensionMismatch(
        "Cannot add a row with $nz indices but $(length(rval)) non-zeros"
    ))

    # Go through through rval to check all coeffs are finite and remove zeros.
    _rind = Vector{Int}(undef, nz)
    _rval = Vector{T}(undef, nz)
    _nz = 0
    for (j, aij) in zip(rind, rval)
        if !iszero(aij)
            isfinite(aij) || error("Invalid row coefficient: $(aij)")
            _nz += 1
            _rind[_nz] = j
            _rval[_nz] = aij
        end
    end
    resize!(_rind, _nz)
    resize!(_rval, _nz)

    # TODO: combine dupplicate indices

    # Increment row counter
    pb.ncon += 1
    push!(pb.lcon, l)
    push!(pb.ucon, u)
    push!(pb.con_names, name)

    # Create new row
    if issorted
        row = Row{T}(_rind, _rval)
    else
        # Sort indices first
        p = sortperm(_rind)
        row = Row{T}(_rind[p], _rval[p])
    end
    push!(pb.arows, row)

    # Update column coefficients
    for (j, aij) in zip(_rind, _rval)
        push!(pb.acols[j].nzind, pb.ncon)
        push!(pb.acols[j].nzval, aij)
    end

    # Done
    return pb.ncon
end

"""
    add_variable!(pb, cind, cval, obj, l, u, [name])

Add one variable to the problem.

# Arguments
* `pb::ProblemData{T}`: the problem to which the new column is added
* `cind::Vector{Int}`: row indices in the new column
* `cval::Vector{T}`: non-zero values in the new column
* `obj::T`: objective coefficient
* `l::T`: column lower bound
* `u::T`: column upper bound
* `name::String`: column name (defaults to `""`)
* `issorted::Bool`: indicates whether the column indices are already issorted.
"""
function add_variable!(pb::ProblemData{T},
    cind::Vector{Int}, cval::Vector{T},
    obj::T, l::T, u::T,
    name::String="";
    issorted::Bool=false
)::Int where{T}
    # Sanity checks
    nz = length(cind)
    nz == length(cval) || throw(DimensionMismatch(
        "Cannot add a column with $nz indices but $(length(cval)) non-zeros"
    ))

    # Go through through cval to check all coeffs are finite and remove zeros.
    _cind = Vector{Int}(undef, nz)
    _cval = Vector{T}(undef, nz)
    _nz = 0
    for (j, aij) in zip(cind, cval)
        if !iszero(aij)
            isfinite(aij) || error("Invalid column coefficient: $(aij)")
            _nz += 1
            _cind[_nz] = j
            _cval[_nz] = aij
        end
    end
    resize!(_cind, _nz)
    resize!(_cval, _nz)

    # Increment column counter
    pb.nvar += 1
    push!(pb.lvar, l)
    push!(pb.uvar, u)
    push!(pb.obj, obj)
    push!(pb.var_names, name)

    # TODO: combine dupplicate indices

    # Create a new column
    if issorted
        col = Col{T}(_cind, _cval)
    else
        # Sort indices
        p = sortperm(_cind)
        col = Col{T}(_cind[p], _cval[p])
    end
    push!(pb.acols, col)

    # Update row coefficients
    for (i, aij) in zip(_cind, _cval)
        push!(pb.arows[i].nzind, pb.nvar)
        push!(pb.arows[i].nzval, aij)
    end

    # Done
    return pb.nvar
end

"""
    load_problem!(pb, )

Load entire problem.
"""
function load_problem!(pb::ProblemData{T},
    name::String,
    objsense::Bool, obj::Vector{T}, obj0::T,
    A::SparseMatrixCSC,
    lcon::Vector{T}, ucon::Vector{T},
    lvar::Vector{T}, uvar::Vector{T},
    con_names::Vector{String}, var_names::Vector{String}
) where{T}
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
    pb.acols = Vector{Col{T}}(undef, nvar)
    pb.arows = Vector{Row{T}}(undef, ncon)
    for j in 1:nvar
        col = A[:, j]
        pb.acols[j] = Col{T}(col.nzind, col.nzval)
    end

    At = sparse(A')
    for i in 1:ncon
        row = At[:, i]
        pb.arows[i] = Row{T}(row.nzind, row.nzval)
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
function delete_constraint!(pb::ProblemData{T}, rind::Int) where{T}
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
function delete_constraints!(pb::ProblemData{T}, rinds) where{T}
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
function delete_variable!(pb::ProblemData{T}, cind::Int) where{T}
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
function delete_variables!(pb::ProblemData{T}, cinds) where{T}
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
* `pb::ProblemData{T}`: the problem whose coefficient
* `i::Int`: row index
* `j::Int`: column index
* `v::T`: coefficient value
"""
function set_coefficient!(pb::ProblemData{T}, i::Int, j::Int, v::T) where{T}
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
    _set_coefficient!(roc::RowOrCol{T}, ind::Int, v::T)

Set coefficient to value `v`.
"""
function _set_coefficient!(roc::RowOrCol{T}, ind::Int, v::T) where{T}
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
