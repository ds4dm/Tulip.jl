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
    lyp::Vector{Tv}
    uyp::Vector{Tv}
    lym::Vector{Tv}
    uym::Vector{Tv}
    lsp::Vector{Tv}
    usp::Vector{Tv}
    lsm::Vector{Tv}
    usm::Vector{Tv}

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
        ps.nzrow = Vector{Int}(undef, ps.nrow)
        ps.nzcol = Vector{Int}(undef, ps.ncol)
        for (j, col) in enumerate(pb.acols)
            ps.nzcol[j] = length(col.nzind)
        end
        for (i, row) in enumerate(pb.arows)
            ps.nzrow[i] = length(row.nzind)
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
        ps.lyp = zeros(Tv, ps.nrow)
        ps.uyp = Vector{Tv}(undef, ps.nrow)
        ps.lym = zeros(Tv, ps.nrow)
        ps.uym = Vector{Tv}(undef, ps.nrow)
        ps.lsp = zeros(Tv, ps.ncol)
        ps.usp = Vector{Tv}(undef, ps.ncol)
        ps.lsm = zeros(Tv, ps.ncol)
        ps.usm = Vector{Tv}(undef, ps.ncol)
        for (i, (lc, uc)) in enumerate(zip(ps.lrow, ps.urow))
            ps.uyp[i] = (lc == Tv(-Inf)) ? zero(Tv) : Tv(Inf)
            ps.uym[i] = (uc == Tv( Inf)) ? zero(Tv) : Tv(Inf)
        end
        for (j, (lv, uv)) in enumerate(zip(ps.lcol, ps.ucol))
            ps.usp[j] = (lv == Tv(-Inf)) ? zero(Tv) : Tv(Inf)
            ps.usm[j] = (uv == Tv( Inf)) ? zero(Tv) : Tv(Inf)
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
function extract_reduced_problem!(dat::PresolveData{Tv}) where{Tv<:Real}
    
    pb = ProblemData{Tv}()

    pb.ncon = sum(dat.rowflag)
    pb.nvar = sum(dat.colflag)

    pb.objsense = dat.objsense
    if pb.objsense
        pb.obj0 = dat.obj0
        pb.obj = dat.obj[dat.colflag]
    else
        pb.obj0 = -dat.obj0
        pb.obj = -dat.obj[dat.colflag]
    end

    # Primal bounds
    pb.lvar = dat.lcol[dat.colflag]
    pb.uvar = dat.ucol[dat.colflag]
    pb.lcon = dat.lrow[dat.rowflag]
    pb.ucon = dat.urow[dat.rowflag]

    # Extract new rows
    pb.arows = Vector{Row{Tv}}(undef, pb.ncon)
    inew = 0
    for (iold, row) in enumerate(dat.pb0.arows)
        dat.rowflag[iold] || continue
        
        inew += 1
        # Compute new row
        rind = Vector{Int}(undef, dat.nzrow[iold])
        rval = Vector{Tv}(undef, dat.nzrow[iold])

        k = 0
        for (jold, aij) in zip(row.nzind, row.nzval)
            dat.colflag[jold] || continue
            # Set new coefficient
            k += 1
            rind[k] = dat.new_var_idx[jold]
            rval[k] = aij
        end

        # Set new row
        pb.arows[inew] = Row{Tv}(rind, rval)
    end

    # Extract new columns
    pb.acols = Vector{Col{Tv}}(undef, pb.nvar)
    jnew = 0
    for (jold, col) in enumerate(dat.pb0.acols)
        dat.colflag[jold] || continue
        
        jnew += 1
        # Compute new row
        cind = Vector{Int}(undef, dat.nzcol[jold])
        cval = Vector{Tv}(undef, dat.nzcol[jold]) 

        k = 0
        for (iold, aij) in zip(col.nzind, col.nzval)
            dat.rowflag[iold] || continue
            # Set new coefficient
            k += 1
            cind[k] = dat.new_con_idx[iold]
            cval[k] = aij
        end

        # Set new column
        pb.acols[jnew] = Col{Tv}(cind, cval)
    end

    # Variable and constraint names
    # TODO: we don't need these
    pb.var_names = dat.pb0.var_names[dat.colflag]
    pb.con_names = dat.pb0.con_names[dat.rowflag]

    # Scaling
    rscale = zeros(Tv, dat.nrow)
    cscale = zeros(Tv, dat.ncol)

    # Compute norm of each row and column
    # TODO: use a parameter p and do norm(.., p)
    for (i, row) in enumerate(pb.arows)
        r = norm(row.nzval, 2)
        rscale[i] = r > zero(Tv) ? r : one(Tv)
    end
    for (j, col) in enumerate(pb.acols)
        r = norm(col.nzval, 2)
        cscale[j] = r > zero(Tv) ? r : one(Tv)
    end

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
    dat.row_scaling = rscale
    dat.col_scaling = cscale

    # Done
    dat.pb_red = pb
    return nothing
end

include("empty_row.jl")
include("empty_column.jl")
include("fixed_variable.jl")
include("row_singleton.jl")


"""
    postsolve!

Perform post-solve.
"""
function postsolve!(sol::Solution{Tv}, sol_::Solution{Tv}, lp::PresolveData{Tv}) where{Tv}

    # Check dimensions
    (sol_.m, sol_.n) == (lp.nrow, lp.ncol) || error(
        "Inner solution has size $((sol_.m, sol_.n)) but presolved problem has size $((lp.nrow, lp.ncol))"
    )
    (sol.m, sol.n) == (lp.pb0.ncon, lp.pb0.nvar) || error(
        "Solution has size $((sol.m, sol.n)) but original problem has size $((lp.pb0.ncon, lp.pb0.nvar))"
    )

    # Copy solution status and objective values
    sol.primal_status = sol_.primal_status
    sol.dual_status = sol_.dual_status
    sol.z_primal = sol_.z_primal
    sol.z_dual = sol_.z_dual

    # Extract and un-scale inner solution components
    # TODO: create a PresolveTransformation for scaling
    for (j_, j) in enumerate(lp.old_var_idx)
        sol.x[j] = sol_.x[j_] / lp.col_scaling[j_]
        sol.s_lower[j] = sol_.s_lower[j_] * lp.col_scaling[j_]
        sol.s_upper[j] = sol_.s_upper[j_] * lp.col_scaling[j_]
    end
    for (i_, i) in enumerate(lp.old_con_idx)
        sol.y_lower[i] = sol_.y_lower[i_] / lp.row_scaling[i_]
        sol.y_upper[i] = sol_.y_upper[i_] / lp.row_scaling[i_]
    end

    # Reverse transformations
    for op in Iterators.reverse(lp.ops)
        postsolve!(sol, sol_, op)
    end

    # Compute row primals
    for (i, row) in enumerate(lp.pb0.arows)
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
function presolve!(lp::PresolveData{Tv}) where{Tv<:Real}
    tstart = time()

    # TODO: check bound consistency on all rows/columns
    st = bounds_consistency_checks!(lp)

    if lp.status == Trm_PrimalInfeasible
        # Problem is infeasible => stop here
        compute_index_mapping!(lp)
        return lp.status
    end

    # I. Remove all fixed variables, empty rows and columns
    # remove_fixed_variables!(lp)
    remove_empty_rows!(lp)
    remove_empty_columns!(lp)

    # TODO: check status for potential early return
    if lp.status != Trm_Unknown
        compute_index_mapping!(lp)
        return lp.status
    end

    # Identify row singletons
    lp.row_singletons = [i for (i, nz) in enumerate(lp.nzrow) if lp.rowflag[i] && nz == 1]
    
    # II. Passes
    lp.updated = true
    npasses = 0  # TODO: maximum number of passes
    while lp.updated && lp.status == Trm_Unknown
        npasses += 1
        lp.updated = false
        @info "Presolve pass $npasses"

        bounds_consistency_checks!(lp)
        lp.status == Trm_Unknown || break

        # Remove all fixed variables
        remove_fixed_variables!(lp)
        lp.status == Trm_Unknown || break

        # Remove row singletons
        remove_row_singletons!(lp)
        lp.status == Trm_Unknown || break

        # Remove forcing & dominated constraints
        # remove_forcing_constraints!(lp)

        # Remove free and implied free column singletons
        # remove_free_column_singletons!(lp)

        # TODO: remove column singleton with doubleton equation

        # Dual reductions
        # remove_dominated_columns!(lp)

        # TODO: fixed variables, empty row/column should be eliminated on the spot
        remove_empty_rows!(lp)
        lp.status == Trm_Unknown || break

        remove_empty_columns!(lp)
        lp.status == Trm_Unknown || break

    end

    @info("Presolved model stats:",
        sum(lp.rowflag), sum(lp.colflag)
    )

    # TODO: check problem dimensions and declare optimality if problem is empty
    if lp.nrow == 0 && lp.ncol == 0
        # Problem is empty: declare optimality now
        lp.status = Trm_Optimal

        # Resize solution
        resize!(lp.solution, 0, 0)
        lp.solution.primal_status = Sln_Optimal
        lp.solution.dual_status = Sln_Optimal
        lp.solution.z_primal = lp.obj0
        lp.solution.z_dual = lp.obj0
    end
    
    # Old <-> new index mapping
    compute_index_mapping!(lp)

    # TODO: extract reduced problem (?)

    # Done.
    @info "Presolve time "*@sprintf("%.2fs", time() - tstart)
    return lp.status
end

function compute_index_mapping!(lp::PresolveData)
    lp.new_con_idx = Vector{Int}(undef, lp.pb0.ncon)
    lp.new_var_idx = Vector{Int}(undef, lp.pb0.nvar)
    lp.old_con_idx = Vector{Int}(undef, lp.nrow)
    lp.old_var_idx = Vector{Int}(undef, lp.ncol)

    inew = 0
    for iold in 1:lp.pb0.ncon
        if lp.rowflag[iold]
            inew += 1
            lp.new_con_idx[iold] = inew
            lp.old_con_idx[inew] = iold
        else
            lp.new_con_idx[iold] = 0
        end
    end
    jnew = 0
    for jold in 1:lp.pb0.nvar
        if lp.colflag[jold]
            jnew += 1
            lp.new_var_idx[jold] = jnew
            lp.old_var_idx[jnew] = jold
        else
            lp.new_var_idx[jold] = 0
        end
    end

    return nothing
end

"""
    bounds_consistency_checks(lp)

Check that all primal & dual bounds are consistent.

TODO: If not, declare primal/dual infeasibility and extract ray.
"""
function bounds_consistency_checks!(lp::PresolveData{Tv}) where{Tv}
    # Check primal bounds
    for (i, (l, u)) in enumerate(zip(lp.lrow, lp.urow))
        if lp.rowflag[i] && l > u
            # Problem is primal infeasible
            @info "Row $i is primal infeasible"
            lp.status = Trm_PrimalInfeasible
            lp.updated = true

            # Farkas ray: y⁺_i = y⁻_i = 1 (any > 0 value works)
            resize!(lp.solution, lp.nrow, lp.ncol)
            lp.solution.x .= zero(Tv)
            lp.solution.y_lower .= zero(Tv)
            lp.solution.y_upper .= zero(Tv)
            lp.solution.s_lower .= zero(Tv)
            lp.solution.s_upper .= zero(Tv)

            lp.solution.primal_status = Sln_Unknown
            lp.solution.dual_status = Sln_InfeasibilityCertificate
            lp.solution.y_lower[i] = one(Tv)
            lp.solution.y_upper[i] = one(Tv)

            return
        end
    end
    for (j, (l, u)) in enumerate(zip(lp.lcol, lp.ucol))
        if lp.colflag[j] && l > u
            # Primal is primal infeasible
            @info "Column $j is primal infeasible"
            lp.status = Trm_PrimalInfeasible
            lp.updated = true

            # Farkas ray
            resize!(lp.solution, lp.nrow, lp.ncol)
            lp.solution.x .= zero(Tv)
            lp.solution.y_lower .= zero(Tv)
            lp.solution.y_upper .= zero(Tv)
            lp.solution.s_lower .= zero(Tv)
            lp.solution.s_upper .= zero(Tv)

            lp.solution.primal_status = Sln_Unknown
            lp.solution.dual_status = Sln_InfeasibilityCertificate
            lp.solution.s_lower[j] = one(Tv)
            lp.solution.s_upper[j] = one(Tv)
            
            return
        end
    end

    # TODO: Check dual bounds

    return nothing
end

"""
    remove_empty_rows!(lp::PresolveData)

Remove all empty rows.

Called once at the beginning of the presolve procedure.
If an empty row is created later, it is removed on the spot.
"""
function remove_empty_rows!(lp::PresolveData{Tv}) where{Tv}
    nempty = 0
    for i in 1:lp.pb0.ncon
        (lp.rowflag[i] && (lp.nzrow[i] == 0)) || continue
        @debug "Remove empty row $i"

        remove_empty_row!(lp, i)
    end
    return nothing
end

"""
    remove_empty_columns!(lp::PresolveData)

Remove all empty columns.

Called once at the beginning of the presolve procedure.
If an empty column is created later, it is removed on the spot.
"""
function remove_empty_columns!(lp::PresolveData{Tv}) where{Tv}
    for j in 1:lp.pb0.nvar
        remove_empty_column!(lp, j)
        lp.status == Trm_Unknown || break
    end
    return nothing
end

"""
    remove_fixed_variables!(lp::PresolveData)

Remove all fixed variables.
"""
function remove_fixed_variables!(lp::PresolveData{Tv}) where{Tv}
    for (j, flag) in enumerate(lp.colflag)
        flag || continue
        remove_fixed_variable!(lp, j)
    end
    return nothing
end

function remove_row_singletons!(lp::PresolveData{Tv}) where{Tv}
    nsingletons = 0
    for i in lp.row_singletons
        remove_row_singleton!(lp, i)
    end
    lp.row_singletons = Int[]
    return nothing
end