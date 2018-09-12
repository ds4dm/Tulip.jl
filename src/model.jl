"""
    Model
    
Data structure for a model.

# Attributes
"""
mutable struct Model
    #=======================================================
        Optimization environment
    =======================================================#

    env::TulipEnv                   # Environment
    sol_status::SolutionStatus      # Solution status
    
    time_total::Float64     # Total time spent by the optimization
    num_bar_iter::Int       # Total number of Interior-Point iterations
    primal_bound::Float64   # Best known primal bound
    dual_bound::Float64     # Best known dual bound

    #=======================================================
        Problem data
    =======================================================#

    num_var::Int              # Number of variables (total)
    n_var_ub::Int           # Number of upper-bounded variables
    num_constr::Int           # Number of constraints

    A::AbstractMatrix{Float64}      # Constraints matrix
    b::AbstractVector{Float64}      # Right-hand side of equality constraints
    c::AbstractVector{Float64}      # Objective vector
    uind::AbstractVector{Int}       # Indices of upper-bounded variables
    uval::AbstractVector{Float64}   # Values of upper bounds


    #=======================================================
        Book-keeping
    =======================================================#
    x::Vector{Float64}      # Vector of primal variables (original variables)
    w::Vector{Float64}      # Vector of primal variables (upper bound slack)
    y::Vector{Float64}      # Vector of dual variables (equality constraints)
    s::Vector{Float64}      # Vector of reduced costs of `x`
    z::Vector{Float64}      # Vector of reduced costs of `w`
    t::Base.RefValue{Float64}    # Artificial homogeneous variable
    k::Base.RefValue{Float64}    # Artificial homogeneous variable
    μ::Base.RefValue{Float64}    # Current barrier parameter

    rp::Vector{Float64}     # Vector of primal residuals `Ax - b`
    rd::Vector{Float64}     # Vector of dual residuals `A'y + s - c`
    ru::Vector{Float64}     # Vector of primal residuals 'x + w - u`
    rg::Base.RefValue{Float64}   # Residual for optimality gap
    rxs::Vector{Float64}    # Right-hand side for complimentarity product `x*s`
    rwz::Vector{Float64}    # Right-hand side for complimentarity product `w*z`
    rtk::Base.RefValue{Float64}  # Right-hand side for complimentarity product `t*k`

    #=======================================================
        Model constructor
    =======================================================#  
    function Model(
        env::TulipEnv,
        A::AbstractMatrix,
        b::AbstractVector,
        c::AbstractVector,
        uind::AbstractVector{Int},
        uval::AbstractVector,
    )
        m = new()

        m.env = env

        # Dimension check
        num_constr = size(A, 1)
        num_var = size(A, 2)
        num_constr == size(b, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but b has size $(size(b))"
        ))
        num_var == size(c, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but c has size $(size(c))"
        ))
        n_var_ub = size(uind, 1)
        n_var_ub == size(uval, 1) || throw(DimensionMismatch(
            "uind has size $(size(uind)) but uval has size $(size(uval))"
        ))
        n_var_ub <= num_var || throw(DimensionMismatch(
            "Too many upper bounds were specified"
        ))
        if n_var_ub > 0
            uind[end] <= num_var || throw(DimensionMismatch(
                "Got upper bound for var $(uind[end])>$(num_var)"
            ))
        end

        # Copy problem data
        m.num_var = num_var
        m.n_var_ub = n_var_ub
        m.num_constr = num_constr
        m.A = copy(A)
        m.b = copy(b)
        m.c = copy(c)
        m.uind = copy(uind)
        m.uval = copy(uval)

        # Allocate memory for optimization algo
        m.x = Vector{Float64}(num_var)
        m.w = Vector{Float64}(n_var_ub)
        m.y = Vector{Float64}(num_constr)
        m.s = Vector{Float64}(num_var)
        m.z = Vector{Float64}(n_var_ub)
        m.t = Ref(1.0)
        m.k = Ref(1.0)
        m.μ = Ref(1.0)

        m.rp = Vector{Float64}(num_constr)
        m.rd = Vector{Float64}(num_var)
        m.ru = Vector{Float64}(n_var_ub)
        m.rg = Ref(1.0)
        m.rxs = Vector{Float64}(num_var)
        m.rwz = Vector{Float64}(n_var_ub)
        m.rtk = Ref(1.0)

        # Book-keeping stuff
        m.sol_status = Unknown

        return m
    end
end


#=======================================================
    Constructors
=======================================================#

"""
    Model(A, b, c, uind, uval)

Construct a model with upper bounds on the specified variables.

    Model(A, b, c)

Construct a model with no upper bounds.

Model()

Construct an empty model.
"""
Model(A, b, c, uind, uval) = Model(TulipEnv(), A, b, c, uind, uval)
Model(A::Ta, b, c, uind, uval) where{Tv<:Real, Ta<:Matrix{Tv}} = Model(sparse(float(A)), b, c, uind, uval)
# Model(env, A::Ta, b, c, uind, uval) where{Tv<:Real, Ta<:Matrix{Tv}} = Model(env, sparse(float(A)), b, c, uind, uval)

Model(A, b, c) = Model(A, b, c, Vector{Int}(0,), Vector{Float64}(0,))

Model() = Model(spzeros(0, 0), Vector{Float64}(0,), Vector{Float64}(0,))

#=======================================================
    Model interface

Retrieve and modify problem-related information
=======================================================#

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

"""
    getobjectivecoeffs(m::Model)

Return the objective coefficients.
"""
getobjectivecoeffs(m::Model) = copy(m.c)

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
    getvarlowerbounds(m::Model)

Return lower bounds on the variables.
"""
getvarlowerbounds(m::Model) = zeros(m.num_var)

"""
    getvarupperbounds(m::Model)

Return upper bounds on the variables. If a given variable has no explicit upper
    bound, the returned value is `Inf`.
"""
function getvarupperbounds(m::Model)
    ub = fill(Inf, m.num_var)
    ub[m.uind] = m.uval
    return ub
end

"""
    setvarupperbounds!(m::Model, ub)

Set upperbounds on the variables
"""
function setvarupperbounds!(m::Model, ub)
    error("Wrong argument type: ub must be a real-valued vector")
    return nothing
end
function setvarupperbounds!(m::Model, ub::SparseVector{Tv, Ti}) where{Tv<:Real, Ti<:Integer}
    # Dimension check
    size(ub, 1) == m.num_var || throw(DimensionMismatch(
        "ub has size $(size(c, 1)) but model has $(m.num_var) variables"
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
    colvals::AbstractVector{Tv},
    l::Real,
    u::Real,
    objcoef::Real
) where Tv<:Real

    # Dimension check
    m.num_constr == size(colvals, 1) || throw(DimensionMismatch(
        "Column has $(size(col, 1)) coeffs but model has $(m.num_constr) constraints"
    ))
    l <= u || error("Upper-bound must be greater than lower bound")

    if l == -Inf && u == Inf
        # free variable: add positive and negative parts
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, colvals)
        m.num_var += 1

        insert!(m.c, k, objcoef)

        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, -colvals)
        m.num_var += 1

        insert!(m.c, k, -objcoef)


    elseif l == -Inf && u < Inf
        # upperbounded, no lower bound
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, -colvals)
        # insert!(m.x, k, 1.0)
        m.num_var += 1

        # Update objective
        insert!(m.c, k, -objval)

        # Update right-hand side
        m.b .+= u * colvals

    elseif -Inf < l && u < Inf
        # lower- and upper-bounded
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, colvals)
        # insert!(m.x, k, 1.0)
        m.num_var += 1

        # Update objective
        insert!(m.c, k, objcoef)

        # Update right-hand side
        m.b .-= l * colvals

        # Update upper bounds
        insert!(m.uind, k, k)
        insert!(m.uind, k, u-l)
        m.n_var_ub += 1

    else
        # lower-bounded
        # Add column
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, colvals)
        m.num_var += 1

        # Update objective
        insert!(m.c, k, objcoef)

        # Update right-hand side
        m.b .-= l * colvals
    end

    m.status = :Built

    return nothing
end
addvar!(m::Model, constridx, constrcoef, l, u, objcoef) = addvar!(m, sparsevec(constridx, constrcoef, m.num_constr), l, u, objcoef)
addvar!(m::Model, col::AbstractVector{Tv}, objcoef::Real) where Tv<:Real = addvar!(m, col, 0.0, Inf, objcoef)

"""
    getconstrlowerbound

Return lower bound on constraint.
"""
getconstrlowerbounds(m::Model) = copy(m.b)

"""
    getconstrupperbound

Return upper bound on constraint.
"""
getconstrupperbounds(m::Model) = copy(m.b)

"""
    addconstr!(m::Model, rowvals, rhs)

Add a constraint to the model.
"""
function addconstr!(m::Model, row::AbstractVector, rhs::Real)

    # Dimension checks
    # length(row) == length(rhs) || throw(DimensionMismatch(
    #     "Adding $(length(row)) constraints but rhs has $(length(rhs)) elements"
    # ))
    m.num_var == length(row) || throw(DimensionMismatch(
        "Row has $(length(row)) coefs but model has $(m.num_var) variables"
    ))
    min, max = extrema(rhs)
    -Inf < min <= max < Inf || error("Right-hand side must have finite value")

    # Add constraint
    m.num_constr += 1
    m.b = vcat(m.b, rhs)
    m.A = vcat(m.A, row')

    return nothing
end

addconstr!(m, varidx, rowval, rhs) = addconstr!(m, sparsevec(varidx, rowval, m.num_var), rhs)



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
    getobjectivevalue(m::Model)

Return objective value of the best known solution
"""
getobjectivevalue(m::Model) = dot(m.c, m.x) / m.t.x

"""
    getdualbound(m::Model)

Return dual bound on the obective value. Returns a lower (resp. upper) bound
if the problem is a minimization (resp. maximization).
"""
getdualbound(m::Model) = (dot(m.b, m.y) - dot(m.uval, m.z)) / m.t.x

"""
    getobjectivedualgap(m::Model)

Return the duality gap. 
"""
getobjectivedualgap(m::Model) = (dot(m.x, m.s) + dot(m.w, m.z)) / (m.t.x)^2

"""
    getsolution(m::Model)

Return best known (primal) solution to the problem.
"""
getsolution(m::Model) = copy(m.x / m.t.x)

"""
    getconstrduals(m::Model)

Return dual variables associated to linear constraints.
"""
getconstrduals(m::Model) = copy(m.y / m.t.x)

"""
    getreducedcosts(m::Model)

Return reduced costs of primal variables.
"""
getreducedcosts(m::Model) = copy(m.s / m.t.x)

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