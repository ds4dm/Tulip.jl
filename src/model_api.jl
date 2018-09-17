#=======================================================
    Model interface
=======================================================#

#=======================================================
    Model instanciation and information queries
=======================================================#

"""
    loadmodel!()

Load data into existing model.
"""
function loadmodel!(
    m::Model,
    num_var::Int,
    num_constr::Int,
    A::AbstractMatrix,
    obj::AbstractVector,
    var_lb::AbstractVector,
    var_ub::AbstractVector,
    constr_lb::AbstractVector,
    constr_ub::AbstractVector
)
    # Dimension checks
    num_var >= 0 || error("Number of variables must be non-negative.")
    num_constr >= 0 || error("Number of constraints must be non-negative.")
    num_var == length(obj) || throw(DimensionMismatch(
        "Wrong number of elements in obj."
    ))
    num_var == length(var_lb) || throw(DimensionMismatch(
        "Wrong number of elements in var_lb."
    ))
    num_var == length(var_ub) || throw(DimensionMismatch(
        "Wrong number of elements in var_ub."
    ))
    num_var == size(A, 2) || throw(DimensionMismatch(
        "Wrong number of columns in A."
    ))
    num_constr == size(A, 1) || throw(DimensionMismatch(
        "Wrong number of rows in A."
    ))
    num_constr == length(constr_lb) || throw(DimensionMismatch(
        "Wrong number of coefficients in constr_lb"
    ))
    num_constr == length(constr_ub) || throw(DimensionMismatch(
        "Wrong number of coefficients in constr_ub"
    ))

    # Import data
    m.num_var = num_var
    m.num_constr = num_constr
    m.A = copy(A)
    m.obj = copy(obj)
    m.var_lb = copy(var_lb)
    m.var_ub = copy(var_ub)
    m.constr_lb = copy(constr_lb)
    m.constr_ub = copy(constr_ub)

    # De-allocate existing working memory
    m.x = Vector{Float64}(undef, 0)
    m.w = Vector{Float64}(undef, 0)
    m.y = Vector{Float64}(undef, 0)
    m.s = Vector{Float64}(undef, 0)
    m.z = Vector{Float64}(undef, 0)
    m.t = Ref(1.0)
    m.k = Ref(1.0)
    m.μ = Ref(1.0)

    m.rp = Vector{Float64}(undef, 0)
    m.ru = Vector{Float64}(undef, 0)
    m.rd = Vector{Float64}(undef, 0)
    m.rg = Ref(Inf)
    m.rxs = Vector{Float64}(undef, 0)
    m.rwz = Vector{Float64}(undef, 0)
    m.rtk = Ref(Inf)
    
    return true
end

"""
    loadmodel!(m, A, rhs, obj)

Load model in standard form.
"""
loadmodel!(m, A, rhs, obj) = loadmodel!(
    m, size(A, 2), size(A, 1),
    A, obj,
    zeros(size(A, 2)), fill(Inf, size(A, 2)),
    rhs, rhs
)

"""
    getnumvar(m::Model)

Return number of variables in the model. This number does not include artificial
slack variables that are used in the formulation (e.g. related to variable
bounds)
"""
getnumvar(m::Model) = m.num_var

"""
    getnumconstr(m::Model)

Return the number of constraints in the model. This number does not include
    explicit bounds on the variables.
"""
getnumconstr(m::Model) = m.num_constr



#=======================================================
    Problem modification
=======================================================#

"""
    addvar!(m::Model, colval, l, u, objcoeff)

Add one variable to the model.

# Arguments
- `m`: Model
- `colvals`: Coefficients of the column
- `l`: lower bound on the variable
- `u`: Upperbound on the variable
- `objcoef`: Objective coefficient
"""
function addvar!(
    m::Model,
    col::AbstractVector,
    lb::Real,
    ub::Real,
    objcoef::Real
)

    # Check dimensions
    m.num_constr == length(col) || throw(DimensionMismatch(
        "Column has $(length(col)) coeffs but model has $(m.num_constr) constraints."
    ))

    # Check numerical values
    lb <= ub || error("Upper-bound must be greater than lower bound.")
    if m.num_constr > 0
        (min, max) = extrema(col)
        (-Inf < min) && (max < Inf) || error("Coef magnitude if too large.")
    end
    -Inf < objcoef < Inf || error("Obj magnitude is too large.")

    # Update model
    m.A, k = Tulip.TLPLinearAlgebra.addcolumn!(m.A, col)
    insert!(m.obj, k, objcoef)
    insert!(m.var_lb, k, lb)
    insert!(m.var_ub, k, ub)
    m.num_var += 1

    return k
end
addvar!(m::Model, constridx, constrcoef, l, u, objcoef) = addvar!(m, sparsevec(constridx, constrcoef, m.num_constr), l, u, objcoef)
addvar!(m::Model, col::AbstractVector, objcoef) = addvar!(m, col, 0.0, Inf, objcoef)

"""
    addconstr!(m::Model, rowvals, rhs)

Add a constraint to the model.
"""
function addconstr!(
    m::Model,
    row::AbstractVector,
    lb::Real,
    ub::Real
)

    # Check dimensions
    m.num_var == length(row) || throw(DimensionMismatch(
        "Row has $(length(row)) coefs but model has $(m.num_var) variables."
    ))

    # Check numerical values
    lb <= ub || error("Upper-bound must be greater than lower bound.")
    if m.num_constr > 0
        (min, max) = extrema(row)
        (-Inf < min) && (max < Inf) || error("Coef magnitude if too large.")
    end

    # Update model
    m.A, k = Tulip.TLPLinearAlgebra.addrow!(m.A, row)
    insert!(m.constr_lb, k, lb)
    insert!(m.constr_ub, k, ub)
    m.num_constr += 1

    return nothing
end

addconstr!(m, varidx, rowval, rhs) = addconstr!(m, sparsevec(varidx, rowval, m.num_var), rhs)


#=======================================================
    Solution information
=======================================================#
"""
    getsolution(m::Model)

Return best known (primal) solution to the problem.
"""
getprimalsolution(m::Model) = copy(m.x / m.t.x)

"""
    getreducedcosts(m::Model)

Return reduced costs of primal variables.
"""
getreducedcosts(m::Model) = copy(m.s / m.t.x)

"""
    getconstrduals(m::Model)

Return dual variables associated to linear constraints.
"""
getdualsolution(m::Model) = copy(m.y / m.t.x)

"""
    getobjectivevalue(m::Model)

Return objective value of the best known solution
"""
function getobjectivevalue(m::Model)
    s = m.sln_status
    if (
        s == Sln_Optimal
        || s == Sln_PrimalFeasible
        || s == Sln_PrimalDualFeasible
    )
        # A feasible solution exists.
        return dot(m.c, m.x) / m.t.x

    elseif s == Sln_DualInfeasible
        # Problem is unbounded.
        return -Inf

    elseif s == Sln_PrimalInfeasible || s == Sln_PrimalDualInfeasible
        # Problem is infeasible.
        return Inf

    else
        # No solution is available.
        error("No solution available.")
    end
end

"""
    getdualbound(m::Model)

Return dual bound on the obective value. Returns a lower (resp. upper) bound
if the problem is a minimization (resp. maximization).
"""
function getdualbound(m::Model)
    s = m.sln_status
    if (
        s == Sln_Optimal
        || s == Sln_DualFeasible
        || s == Sln_PrimalDualFeasible
    )
        # A feasible solution exists.
        return (dot(m.b, m.y) - dot(m.uval, m.z)) / m.t.x

    elseif s == Sln_DualInfeasible
        # Dual is infeasible.
        return -Inf

    elseif s == Sln_PrimalInfeasible || s == Sln_PrimalDualInfeasible
        # Primal is unbounded.
        return Inf

    else
        # No solution is available.
        error("No dual solution available.")
    end
end

"""
    getobjectivedualgap(m::Model)

Return the duality gap. 
"""
function getobjectivedualgap(m::Model) 
    s = m.sln_status
    if s == Sln_Optimal || s == Sln_PrimalDualFeasible
        # Problem is primal-dual feasible.
        return (dot(m.x, m.s) + dot(m.w, m.z)) / (m.t.x)^2
    else
        return Inf
    end

end


#=======================================================
    Parameter management
=======================================================#


#=======================================================
    Others
=======================================================#
"""
    optimize!(model)

Solve the problem.
"""
function optimize!(::Model) end

"""
    getobjectivecoeffs(m::Model)

Return the objective coefficients.
"""
getobjectivecoeffs(m::Model) = copy(m.obj)

"""
    setobjectivecoeffs!(m::Model, c)

Set new objective coefficients.
"""
function setobjectivecoeffs!(m::Model, c::AbstractVector{T}) where T<:Real

    # Dimension check
    size(c, 1) == m.num_var || throw(DimensionMismatch(
        "c has $(size(c, 1)) coeffs but model has $(m.num_var) variables"
    ))

    m.c = copy(c)
    return nothing
end

"""
    getvarlowerbound(m, i)

Return lower bound on variable `x[i]`.
"""
getvarlowerbound(m::Model, i::Int) = m.var_lb[i]

"""
    getvarupperbound(m, i)

Return upper bound on variable `x[i]`.
"""
getvarupperbound(m::Model, i::Int) = m.var_ub[i]

"""
    getvarlowerbounds(m)

Return lower bounds on all variables.

    getvarlowerbounds(m, idx)

Return lower bounds on specified variables.

    getvarlowerbounds!(m, l)

Overwrite `l` with the lower bounds on variables.

    getvarlowerbounds!(m, idx, l)

Overwrite `l` with the lower bounds on specified variables.
"""
getvarlowerbounds(m::Model) = copy(m.var_lb)

getvarlowerbounds(m::Model, idx::AbstractVector{Int}) = copy(m.var_lb[idx])

getvarlowerbounds!(m::Model, l::Vector) = copy!(l, m.var_lb)

getvarlowerbounds!(
    m::Model,
    l::Vector,
    idx::AbstractVector{Int}
) = copy!(l, m.var_lb[idx])

"""
    getvarupperbounds(m::Model)

Return upper bounds on the variables. If a given variable has no explicit upper
    bound, the returned value is `Inf`.
"""
getvarupperbounds(m::Model) = copy(m.var_ub)
getvarupperbounds!(m::Model, u::Vector) = copy!(u, m.var_ub)

"""
    setvarupperbounds!(m::Model, ub)

Set upperbounds on the variables
"""
function setvarupperbounds!(m::Model, ub::SparseVector{Tv, Ti}) where{Tv<:Real, Ti<:Integer}
    # Dimension check
    size(ub, 1) == m.num_var || throw(DimensionMismatch(
        "ub has size $(size(ub, 1)) but model has $(m.num_var) variables"
    ))
    minimum(ub) >= zero(Tv) || error("Upper bounds must be non-negative")

    m.uind = ub.nzind
    m.uval = ub.nzval
    m.n_var_ub = nnz(ub)
    return nothing
end

function setvarupperbounds!(m::Model, ub::AbstractArray{Tv}) where{Tv<:Real}
    # Dimension check
    size(ub, 1) == m.num_var || throw(DimensionMismatch(
        "ub has size $(size(c, 1)) but model has $(m.num_var) variables"
    ))
    minimum(ub) >= zero(Tv) || error("Upper bounds must be non-negative")

    u_ = ub .< Inf
    m.uind = collect(1:m.num_var)[u_]
    m.uval = ub[u_]
    m.n_var_ub = sum(u_)

    return nothing
end


"""
    getconstrlowerbound(m, j)

Return lower bound on constraint `j`.

    getconstrlowerbounds(m, J)

Return lower bounds on constraints `j` for `j` in `J`.

    getconstrlowerbounds(m)

Return lower bounds on all constraints.
"""
getconstrlowerbound(m::Model, j) = m.constr_lb[j]
getconstrlowerbounds(m::Model, J::AbstractVector{Integer}) = copy(m.constr_lb[J])
getconstrlowerbounds(m::Model) = copy(m.constr_lb)


"""
    getconstrupperbound(m, j)

Return upper bound on constraint `j`.

    getconstrupperbounds(m, J)

Return upper bounds on constraints `j` for `j` in `J`.

    getconstrupperbounds(m)

Return upper bounds on all constraints.
"""
getconstrupperbound(m::Model, j) = m.constr_ub[j]
getconstrupperbounds(m::Model, J::AbstractVector{Integer}) = copy(m.constr_ub[J])
getconstrupperbounds(m::Model) = copy(m.constr_ub)




"""
    getlinearconstrcoeffs(m::Model)

Return the matrix of linear constraints.
"""
getlinearconstrcoeffs(m::Model) = copy(m.A)


#=======================================================
    Solution interface

Retrieve solution-related information
=======================================================#

"""
    getnumbarrieriter(m::Model)

Return number of barrier iterations.
"""
getnumbarrieriter(m::Model) = m.num_bar_iter

"""
    getsolutiontime(m::Model)

Return runtime of the optimizer.
"""
getsolutiontime(m::Model) = m.time_total

"""
    getinfeasibilityray(m::Model)

Retrieve infeasibility ray when problem is proven infeasible.
"""
getinfeasibilityray(m::Model) = copy(m.y ./ m.t.x)

"""
    getunboundedray(m::Model)

Retrieve unbounded ray when model is proven unbounded.
"""
getunboundedray(m::Model) = copy(m.x ./ m.t.x)