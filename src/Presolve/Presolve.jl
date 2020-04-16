"""
    PresolveTransformation{Tv}

Abstract type for pre-solve transformations.
"""
abstract type PresolveTransformation{Tv} end

"""
    PresolveData{Tv}

Stores information about an LP in the form
```
    min     c'x + c0
    s.t.    lr ⩽ Ax ⩽ ur
            lc ⩽  x ⩽ uc
```
whose dual writes
```
    max     lr'y⁺ - ur'y⁻ + lc's⁺ - uc's⁻
    s.t.     A'y⁺ -  A'y⁻ +    s⁺ -    s⁻ = c
               y⁺,     y⁻,     s⁺,     s⁻ ⩾ 0
```
"""
mutable struct PresolveData{Tv<:Real}
    updated::Bool
    status::TerminationStatus

    # Original problem
    pb0::ProblemData{Tv}
    # Reduced problem
    # Nothing until the reduced problem is extracted
    pb_red::Union{Nothing, ProblemData{Tv}}
    solution::Solution{Tv}  # only used if presolve solves the problem

    # Presolved data

    # Active rows and columns
    rowflag::Vector{Bool}
    colflag::Vector{Bool}

    # Non-zeros in rows and columns
    nzrow::Vector{Int}
    nzcol::Vector{Int}

    # Objective
    objsense::Bool
    obj::Vector{Tv}
    obj0::Tv

    # Current number of constraints/variables in presolved problem
    nrow::Int
    ncol::Int

    # Primal bounds
    lrow::Vector{Tv}
    urow::Vector{Tv}
    lcol::Vector{Tv}
    ucol::Vector{Tv}

    # Dual bounds
    ly::Vector{Tv}
    uy::Vector{Tv}
    ls::Vector{Tv}
    us::Vector{Tv}

    # Scaling
    row_scaling::Vector{Tv}
    col_scaling::Vector{Tv}
    # TODO: objective and RHS scaling

    # Old <-> new index mapping
    # Instantiated only after pre-solve is performed
    new_con_idx::Vector{Int}
    new_var_idx::Vector{Int}
    old_con_idx::Vector{Int}
    old_var_idx::Vector{Int}

    # Singletons
    row_singletons::Vector{Int}  # Row singletons
    free_col_singletons::Vector{Int}  # (implied) free column singletons

    # TODO: set of transformations for pre-post crush
    ops::Vector{PresolveTransformation{Tv}}

    function PresolveData(pb::ProblemData{Tv}) where{Tv}
        ps = new{Tv}()

        ps.updated = false
        ps.status = Trm_Unknown

        ps.pb0 = pb
        ps.pb_red = nothing
        ps.solution = Solution{Tv}(pb.ncon, pb.nvar)

        ps.nrow = pb.ncon
        ps.ncol = pb.nvar

        # All rows and columns are active
        ps.rowflag = trues(ps.nrow)
        ps.colflag = trues(ps.ncol)

        # Number of non-zeros in rows/columns
        ps.nzrow = zeros(Int, ps.nrow)
        ps.nzcol = zeros(Int, ps.ncol)
        for (j, col) in enumerate(pb.acols)
            for (i, aij) in zip(col.nzind, col.nzval)
                ps.nzcol[j] += !iszero(aij)
                ps.nzrow[i] += !iszero(aij)
            end
        end

        # Objective
        ps.objsense = pb.objsense
        if pb.objsense
            ps.obj  = copy(pb.obj)
            ps.obj0 = pb.obj0
        else
            # Maximization problem: negate the objective for pre-solve
            # This will be undone when extracting the reduced problem
            ps.obj  = -copy(pb.obj)
            ps.obj0 = -pb.obj0
        end

        # Copy primal bounds
        ps.lrow = copy(pb.lcon)
        ps.urow = copy(pb.ucon)
        ps.lcol = copy(pb.lvar)
        ps.ucol = copy(pb.uvar)

        # Set dual bounds
        ps.ly = Vector{Tv}(undef, ps.nrow)
        ps.uy = Vector{Tv}(undef, ps.nrow)
        ps.ls = Vector{Tv}(undef, ps.ncol)
        ps.us = Vector{Tv}(undef, ps.ncol)
        for (i, (lc, uc)) in enumerate(zip(ps.lrow, ps.urow))
            ps.ly[i] = (uc == Tv( Inf)) ? zero(Tv) : Tv(-Inf)
            ps.uy[i] = (lc == Tv(-Inf)) ? zero(Tv) : Tv( Inf)
        end
        for (j, (lv, uv)) in enumerate(zip(ps.lcol, ps.ucol))
            ps.ls[j] = (uv == Tv( Inf)) ? zero(Tv) : Tv(-Inf)
            ps.us[j] = (lv == Tv(-Inf)) ? zero(Tv) : Tv( Inf)
        end

        # Scalings
        ps.row_scaling = ones(Tv, ps.nrow)
        ps.col_scaling = ones(Tv, ps.ncol)

        # Index mappings
        ps.new_con_idx = Int[]
        ps.new_var_idx = Int[]
        ps.old_con_idx = Int[]
        ps.old_var_idx = Int[]

        # Singletons
        ps.row_singletons = Int[]
        ps.free_col_singletons = Int[]

        ps.ops = PresolveTransformation{Tv}[]

        return ps
    end
end

# Extract pre-solved problem data, to be passed to the IPM solver
function extract_reduced_problem!(ps::PresolveData{Tv}) where{Tv<:Real}
    
    pb = ProblemData{Tv}()

    pb.ncon = sum(ps.rowflag)
    pb.nvar = sum(ps.colflag)

    pb.objsense = ps.objsense
    if pb.objsense
        pb.obj0 = ps.obj0
        pb.obj = ps.obj[ps.colflag]
    else
        pb.obj0 = -ps.obj0
        pb.obj = -ps.obj[ps.colflag]
    end

    # Primal bounds
    pb.lvar = ps.lcol[ps.colflag]
    pb.uvar = ps.ucol[ps.colflag]
    pb.lcon = ps.lrow[ps.rowflag]
    pb.ucon = ps.urow[ps.rowflag]

    # Extract new rows
    pb.arows = Vector{Row{Tv}}(undef, pb.ncon)
    inew = 0
    for (iold, row) in enumerate(ps.pb0.arows)
        ps.rowflag[iold] || continue
        
        inew += 1
        # Compute new row
        rind = Vector{Int}(undef, ps.nzrow[iold])
        rval = Vector{Tv}(undef, ps.nzrow[iold])

        k = 0
        for (jold, aij) in zip(row.nzind, row.nzval)
            ps.colflag[jold] || continue
            # Set new coefficient
            k += 1
            rind[k] = ps.new_var_idx[jold]
            rval[k] = aij
        end

        # Set new row
        pb.arows[inew] = Row{Tv}(rind, rval)
    end

    # Extract new columns
    pb.acols = Vector{Col{Tv}}(undef, pb.nvar)
    jnew = 0
    for (jold, col) in enumerate(ps.pb0.acols)
        ps.colflag[jold] || continue
        
        jnew += 1
        # Compute new row
        cind = Vector{Int}(undef, ps.nzcol[jold])
        cval = Vector{Tv}(undef, ps.nzcol[jold]) 

        k = 0
        for (iold, aij) in zip(col.nzind, col.nzval)
            ps.rowflag[iold] || continue
            # Set new coefficient
            k += 1
            cind[k] = ps.new_con_idx[iold]
            cval[k] = aij
        end

        # Set new column
        pb.acols[jnew] = Col{Tv}(cind, cval)
    end

    # Variable and constraint names
    # TODO: we don't need these
    pb.var_names = ps.pb0.var_names[ps.colflag]
    pb.con_names = ps.pb0.con_names[ps.rowflag]

    # Scaling
    rscale = zeros(Tv, ps.nrow)
    cscale = zeros(Tv, ps.ncol)

    # Compute norm of each row and column
    # TODO: use a parameter p and do norm(.., p)
    p = 2
    for (i, row) in enumerate(pb.arows)
        r = norm(row.nzval, p)
        rscale[i] = r > zero(Tv) ? r : one(Tv)
    end
    for (j, col) in enumerate(pb.acols)
        r = norm(col.nzval, p)
        cscale[j] = r > zero(Tv) ? r : one(Tv)
    end

    map!(sqrt, cscale, cscale)
    map!(sqrt, rscale, rscale)

    # Rows
    for (i, row) in enumerate(pb.arows)
        # Scale row coefficients
        for (k, j) in enumerate(row.nzind)
            row.nzval[k] /= (rscale[i] * cscale[j])
        end
        # Scale row bounds
        pb.lcon[i] /= rscale[i]
        pb.ucon[i] /= rscale[i]
    end
    # Columns
    for (j, col) in enumerate(pb.acols)
        # Scale column coefficients
        for (k, i) in enumerate(col.nzind)
            col.nzval[k] /= (rscale[i] * cscale[j])
        end
        # Scale objective and variable bounds
        pb.obj[j]  /= cscale[j]
        pb.lvar[j] *= cscale[j]
        pb.uvar[j] *= cscale[j]
    end

    # Record scaling
    @debug "Scaling info" extrema(rscale) extrema(cscale)
    ps.row_scaling = rscale
    ps.col_scaling = cscale

    # Done
    ps.pb_red = pb
    return nothing
end

include("empty_row.jl")
include("empty_column.jl")
include("fixed_variable.jl")
include("row_singleton.jl")
include("forcing_row.jl")
include("free_column_singleton.jl")
include("dominated_column.jl")


"""
    postsolve!

Perform post-solve.
"""
function postsolve!(sol::Solution{Tv}, sol_::Solution{Tv}, ps::PresolveData{Tv}) where{Tv}

    # Check dimensions
    (sol_.m, sol_.n) == (ps.nrow, ps.ncol) || error(
        "Inner solution has size $((sol_.m, sol_.n)) but presolved problem has size $((ps.nrow, ps.ncol))"
    )
    (sol.m, sol.n) == (ps.pb0.ncon, ps.pb0.nvar) || error(
        "Solution has size $((sol.m, sol.n)) but original problem has size $((ps.pb0.ncon, ps.pb0.nvar))"
    )

    # Copy solution status and objective values
    sol.primal_status = sol_.primal_status
    sol.dual_status = sol_.dual_status
    sol.is_primal_ray = sol_.is_primal_ray
    sol.is_dual_ray = sol_.is_dual_ray
    sol.z_primal = sol_.z_primal
    sol.z_dual = sol_.z_dual

    # Extract and un-scale inner solution components
    # TODO: create a PresolveTransformation for scaling
    for (j_, j) in enumerate(ps.old_var_idx)
        sol.x[j] = sol_.x[j_] / ps.col_scaling[j_]
        sol.s_lower[j] = sol_.s_lower[j_] * ps.col_scaling[j_]
        sol.s_upper[j] = sol_.s_upper[j_] * ps.col_scaling[j_]
    end
    for (i_, i) in enumerate(ps.old_con_idx)
        sol.y_lower[i] = sol_.y_lower[i_] / ps.row_scaling[i_]
        sol.y_upper[i] = sol_.y_upper[i_] / ps.row_scaling[i_]
    end

    # Reverse transformations
    for op in Iterators.reverse(ps.ops)
        postsolve!(sol, op)
    end

    # Compute row primals
    for (i, row) in enumerate(ps.pb0.arows)
        sol.Ax[i] = zero(Tv)
        for (j, aij) in zip(row.nzind, row.nzval)
            sol.Ax[i] += aij * sol.x[j]
        end
    end

    # Done
    return nothing
end


"""
    presolve(pb::ProblemData)

Perform pre-solve.
"""
function presolve!(ps::PresolveData{Tv}) where{Tv<:Real}

    # Check bound consistency on all rows/columns
    st = bounds_consistency_checks!(ps)
    ps.status == Trm_PrimalInfeasible && return ps.status

    # I. Remove all fixed variables, empty rows and columns
    # remove_fixed_variables!(ps)
    remove_empty_rows!(ps)
    remove_empty_columns!(ps)

    # TODO: check status for potential early return
    ps.status == Trm_Unknown || return ps.status

    # Identify row singletons
    ps.row_singletons = [i for (i, nz) in enumerate(ps.nzrow) if ps.rowflag[i] && nz == 1]
    
    # II. Passes
    ps.updated = true
    npasses = 0  # TODO: maximum number of passes
    while ps.updated && ps.status == Trm_Unknown
        npasses += 1
        ps.updated = false
        @debug "Presolve pass $npasses" ps.nrow ps.ncol

        bounds_consistency_checks!(ps)
        ps.status == Trm_Unknown || return ps.status
        remove_empty_columns!(ps)
        ps.status == Trm_Unknown || return ps.status
        

        # Remove all fixed variables
        # TODO: remove empty variables as well
        remove_row_singletons!(ps)
        ps.status == Trm_Unknown || return ps.status
        remove_fixed_variables!(ps)
        ps.status == Trm_Unknown || return ps.status

        # Remove forcing & dominated constraints
        remove_row_singletons!(ps)
        ps.status == Trm_Unknown || return ps.status
        remove_forcing_rows!(ps)
        ps.status == Trm_Unknown || return ps.status

        # Remove free and implied free column singletons
        remove_row_singletons!(ps)
        ps.status == Trm_Unknown || return ps.status
        remove_free_column_singletons!(ps)
        ps.status == Trm_Unknown || return ps.status

        # TODO: remove column singleton with doubleton equation

        # Dual reductions
        remove_row_singletons!(ps)
        ps.status == Trm_Unknown || return ps.status
        remove_dominated_columns!(ps)
        ps.status == Trm_Unknown || return ps.status
    end

    remove_empty_columns!(ps)

    @debug("Presolved problem info",
        ps.pb0.ncon, ps.nrow,
        ps.pb0.nvar, ps.ncol,
        sum(ps.nzcol[ps.colflag]), sum(ps.nzrow[ps.rowflag])
    )

    # TODO: check problem dimensions and declare optimality if problem is empty
    if ps.nrow == 0 && ps.ncol == 0
        # Problem is empty: declare optimality now
        ps.status = Trm_Optimal

        # Resize solution
        resize!(ps.solution, 0, 0)
        ps.solution.primal_status = Sln_Optimal
        ps.solution.dual_status = Sln_Optimal
        ps.solution.is_primal_ray = false
        ps.solution.is_dual_ray = false
        ps.solution.z_primal = ps.obj0
        ps.solution.z_dual = ps.obj0
    end
    
    # Old <-> new index mapping
    compute_index_mapping!(ps)

    # TODO: extract reduced problem (?)

    # Done.
    return ps.status
end

function compute_index_mapping!(ps::PresolveData)
    ps.new_con_idx = Vector{Int}(undef, ps.pb0.ncon)
    ps.new_var_idx = Vector{Int}(undef, ps.pb0.nvar)
    ps.old_con_idx = Vector{Int}(undef, ps.nrow)
    ps.old_var_idx = Vector{Int}(undef, ps.ncol)

    inew = 0
    for iold in 1:ps.pb0.ncon
        if ps.rowflag[iold]
            inew += 1
            ps.new_con_idx[iold] = inew
            ps.old_con_idx[inew] = iold
        else
            ps.new_con_idx[iold] = 0
        end
    end
    jnew = 0
    for jold in 1:ps.pb0.nvar
        if ps.colflag[jold]
            jnew += 1
            ps.new_var_idx[jold] = jnew
            ps.old_var_idx[jnew] = jold
        else
            ps.new_var_idx[jold] = 0
        end
    end

    return nothing
end

"""
    bounds_consistency_checks(ps)

Check that all primal & dual bounds are consistent.

TODO: If not, declare primal/dual infeasibility and extract ray.
"""
function bounds_consistency_checks!(ps::PresolveData{Tv}) where{Tv}
    # Check primal bounds
    for (i, (l, u)) in enumerate(zip(ps.lrow, ps.urow))
        if ps.rowflag[i] && l > u
            # Problem is primal infeasible
            @debug "Row $i is primal infeasible"
            ps.status = Trm_PrimalInfeasible
            ps.updated = true

            # Resize problem            
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(Tv)
            ps.solution.y_lower .= zero(Tv)
            ps.solution.y_upper .= zero(Tv)
            ps.solution.s_lower .= zero(Tv)
            ps.solution.s_upper .= zero(Tv)

            # Farkas ray: y⁺_i = y⁻_i = 1 (any > 0 value works)
            ps.solution.primal_status = Sln_Unknown
            ps.solution.dual_status = Sln_InfeasibilityCertificate
            ps.solution.is_primal_ray = false
            ps.solution.is_dual_ray = true
            ps.solution.z_primal = ps.solution.z_dual = Tv(Inf)
            i_ = ps.new_con_idx[i]
            ps.solution.y_lower[i_] = one(Tv)
            ps.solution.y_upper[i_] = one(Tv)
            
            return
        end
    end
    for (j, (l, u)) in enumerate(zip(ps.lcol, ps.ucol))
        if ps.colflag[j] && l > u
            # Primal is primal infeasible
            @debug "Column $j is primal infeasible"
            ps.status = Trm_PrimalInfeasible
            ps.updated = true

            # Resize problem            
            compute_index_mapping!(ps)
            resize!(ps.solution, ps.nrow, ps.ncol)
            ps.solution.x .= zero(Tv)
            ps.solution.y_lower .= zero(Tv)
            ps.solution.y_upper .= zero(Tv)
            ps.solution.s_lower .= zero(Tv)
            ps.solution.s_upper .= zero(Tv)

            # Farkas ray: y⁺_i = y⁻_i = 1 (any > 0 value works)
            ps.solution.primal_status = Sln_Unknown
            ps.solution.dual_status = Sln_InfeasibilityCertificate
            ps.solution.is_primal_ray = false
            ps.solution.is_dual_ray = true
            ps.solution.z_primal = ps.solution.z_dual = Tv(Inf)
            j_ = ps.new_var_idx[j]
            ps.solution.s_lower[j_] = one(Tv)
            ps.solution.s_upper[j_] = one(Tv)
            
            return
        end
    end

    # TODO: Check dual bounds

    return nothing
end

"""
    remove_empty_rows!(ps::PresolveData)

Remove all empty rows.

Called once at the beginning of the presolve procedure.
If an empty row is created later, it is removed on the spot.
"""
function remove_empty_rows!(ps::PresolveData{Tv}) where{Tv}
    nempty = 0
    for i in 1:ps.pb0.ncon
        (ps.rowflag[i] && (ps.nzrow[i] == 0)) || continue
        @debug "Remove empty row $i"

        remove_empty_row!(ps, i)
    end
    return nothing
end

"""
    remove_empty_columns!(ps::PresolveData)

Remove all empty columns.

Called once at the beginning of the presolve procedure.
If an empty column is created later, it is removed on the spot.
"""
function remove_empty_columns!(ps::PresolveData{Tv}) where{Tv}
    for j in 1:ps.pb0.nvar
        remove_empty_column!(ps, j)
        ps.status == Trm_Unknown || break
    end
    return nothing
end

"""
    remove_fixed_variables!(ps::PresolveData)

Remove all fixed variables.
"""
function remove_fixed_variables!(ps::PresolveData{Tv}) where{Tv}
    for (j, flag) in enumerate(ps.colflag)
        flag || continue
        remove_fixed_variable!(ps, j)
    end
    return nothing
end

function remove_row_singletons!(ps::PresolveData{Tv}) where{Tv}
    nsingletons = 0
    for i in ps.row_singletons
        remove_row_singleton!(ps, i)
    end
    ps.row_singletons = Int[]
    return nothing
end

"""
    remove_forcing_rows!

Remove forcing and dominated row
"""
function remove_forcing_rows!(ps::PresolveData)
    for (i, flag) in enumerate(ps.rowflag)
        flag && remove_forcing_row!(ps, i)
    end
    return nothing
end

"""
    remove_free_column_singletons!(ps)

"""
function remove_free_column_singletons!(ps::PresolveData)
    for (j, flag) in enumerate(ps.colflag)
        remove_free_column_singleton!(ps, j)
    end
    return nothing
end

function remove_dominated_columns!(ps::PresolveData{Tv}) where{Tv}
    # Strengthen dual bounds with column singletons
    for (j, (l, u)) in enumerate(zip(ps.lcol, ps.ucol))
        (ps.colflag[j] && ps.nzcol[j] == 1) || continue

        col = ps.pb0.acols[j]
        # Find non-zero index
        nz = 0
        i, aij = 0, zero(Tv)
        for (i_, a_) in zip(col.nzind, col.nzval)
            if ps.rowflag[i_] && !iszero(a_)
                nz += 1; nz <= 1 || break
                i = i_
                aij = a_
            end
        end

        nz == 1 || (@error "Expected singleton but column $j has $nz non-zeros"; continue)
        iszero(aij) && continue  # empty column

        # Strengthen dual bounds
        #=

        =#
        cj = ps.obj[j]
        y_ = cj / aij
        if !isfinite(l) && !isfinite(u)
            # Free variable. Should not happen but handle anyway
            # TODO
            @warn "TODO: dual bounds strengthening for free column singletons"
        elseif isfinite(l) && !isfinite(u)
            # Lower-bounded variable: `aij * yi ≤ cj`
            if aij > zero(Tv)
                # yi ≤ cj / aij
                @debug "Col $j forces y$i <= $y_"
                ps.uy[i] = min(ps.uy[i],  y_)
            else
                # yi ≥ cj / aij
                @debug "Col $j forces y$i >= $y_"
                ps.ly[i] = max(ps.ly[i],  y_)
            end
        
        elseif !isfinite(l) && isfinite(u)
            # Upper-bounded variable: `aij * yi ≥ cj`
            if aij > zero(Tv)
                # yi ≥ cj / aij
                @debug "Col $j forces y$i >= $y_"
                ps.ly[i] = max(ps.ly[i],  y_)
            else
                # yi ≤ cj / aij
                @debug "Col $j forces y$i <= $y_"
                ps.uy[i] = min(ps.uy[i],  y_)
            end
        end

        # TODO: dual feasibility check (?)
    end

    for (j, flag) in enumerate(ps.colflag)
        remove_dominated_column!(ps, j)
        ps.status == Trm_Unknown || break
    end
    return nothing
end