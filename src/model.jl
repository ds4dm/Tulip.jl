"""
    Model
    Data structure for a model

# Attributes
-`n_var::Int`: Number of variables
-`n_con::Int`: Number of constraints
-`A::AbstractMatrix{T1<:Real}`: Constraint matrix
-`b::AbstractVector{T2<:Real}`: Right-hand side of the equality constraints
-`c::AbstractVector{T3<:Real}`: Objective coefficient
-`uind::AbstractVector{Ti<:Integer}`: Indices of upper-bounded variables
-`uval::AbstractVector{T4<:Real}`: Upper bounds on the variables. Only finite
    upper bounds are stored.
-`sol::PrimalDualPoint`: Current solution to the problem (may be infeasible at
    the beginning of the optimization)
-`status::Symbol`: Optimization status
"""
mutable struct Model{Ta<:AbstractMatrix{Float64}}
    #=======================================================
        Optimization environment
    =======================================================#

    env::TulipEnv           # Environment
    status::Symbol          # Optimization status
    numbarrieriter::Int     # Number of barrier iterations
    runtime::Float64      # Elapsed solution time, in seconds


    #=======================================================
        Problem data
    =======================================================#

    n_var::Int              # Number of variables (total)
    n_var_ub::Int           # Number of upper-bounded variables
    n_con::Int              # Number of constraints

    A::Ta                   # Constraints matrix
    b::AbstractVector{Float64}   # Right-hand side of equality constraints
    c::AbstractVector{Float64}   # Objective vector
    uind::AbstractVector{Int}  # Indices of upper-bounded variables
    uval::AbstractVector{Float64}   # Values of upper bounds


    #=======================================================
        Book-keeping
    =======================================================#

    F::Factorization{Float64}    # Factorization object

    x::AbstractVector{Float64}   # Vector of original primal variables
    w::AbstractVector{Float64}   # Vector of primal upper bound slack variables
    y::AbstractVector{Float64}   # Vector of dual variables for equality constraints
    s::AbstractVector{Float64}   # Vector of reduced costs of `x`
    z::AbstractVector{Float64}   # Vector of reduced costs of `w`

    rb::AbstractVector{Float64}  # Vector of primal residuals `Ax - b`
    rc::AbstractVector{Float64}  # Vector of dual residuals `A'y + s - c`
    ru::AbstractVector{Float64}  # Vector of primal residuals 'x + w - u``
    rxs::AbstractVector{Float64} # Right-hand side for omplimentarity product `x*s`
    rwz::AbstractVector{Float64} # Right-hand side for omplimentarity product `w*z`


    #=======================================================
        Model constructor
    =======================================================#  
    function Model(
        env::TulipEnv,
        A::Ta,
        b::AbstractVector,
        c::AbstractVector,
        uind::AbstractVector{Int},
        uval::AbstractVector,
        x::AbstractVector,
        w::AbstractVector,
        y::AbstractVector,
        s::AbstractVector,
        z::AbstractVector,
        rb::AbstractVector,
        rc::AbstractVector,
        ru::AbstractVector,
        rxs::AbstractVector,
        rwz::AbstractVector
    ) where{Ta<:AbstractMatrix{Float64}}
        m = new{Ta}()

        m.env = env

        # Dimension check
        n_con = size(A, 1)
        n_var = size(A, 2)
        n_con == size(b, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but b has size $(size(b))"
        ))
        n_var == size(c, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but c has size $(size(c))"
        ))
        n_var_ub = size(uind, 1)
        n_var_ub == size(uval, 1) || throw(DimensionMismatch(
            "uind has size $(size(uind)) but uval has size $(size(uval))"
        ))
        n_var_ub <= n_var || throw(DimensionMismatch(
            "Too many upper bounds were specified"
        ))
        if n_var_ub > 0
            uind[end] <= n_var || throw(DimensionMismatch(
                "Got upper bound for var $(uind[end])>$(n_var)"
            ))
        end

        m.n_var = n_var
        m.n_var_ub = n_var_ub
        m.n_con = n_con
        m.A = copy(A)
        m.b = copy(b)
        m.c = copy(c)
        m.uind = copy(uind)
        m.uval = copy(uval)

        m.x = copy(x)
        m.w = copy(w)
        m.y = copy(y)
        m.s = copy(s)
        m.z = copy(z)

        m.rb = copy(rb)
        m.rc = copy(rc)
        m.ru = copy(ru)
        m.rxs = copy(rxs)
        m.rwz = copy(rwz)

        m.status = :Built
        m.numbarrieriter = 0
        m.runtime = 0.0

        return m
    end
end

#=======================================================
    Constructors
=======================================================#

"""
    Model(A, b, c, uind, uval)

Construct a model with upper bounds on the specified variables.

    Model()

Construct an empty model.

    Model(A, b, c)

Construct a model with no upper bounds.
"""
function Model(
    env::TulipEnv,
    A::Ta,
    b::AbstractVector{T2},
    c::AbstractVector{T3},
    uind::AbstractVector{Ti},
    uval::AbstractVector{T4}
) where{T1<:Real, Ta<:AbstractMatrix{T1}, T2<:Real, T3<:Real, T4<:Real, Ti<:Integer}
    
    (m, n) = size(A)
    p = size(uind, 1)

    model = Model(
        env,
        float(A), b, c, uind, uval,
        # initial solution
        ones(n),  # x
        ones(p),  # w
        zeros(m), # y
        ones(n),  # s
        ones(p),  # w
        # residuals
        fill(Inf, m),  # rb
        fill(Inf, n),  # rc
        fill(Inf, p),  # ru
        ones(n),  # rxs
        ones(p)   # rwz
    )
    return model
end


Model(A, b, c, uind, uval) = Model(TulipEnv(), A, b, c, uind, uval)
Model(A::Ta, b, c, uind, uval) where{Tv<:Real, Ta<:Matrix{Tv}} = Model(sparse(float(A)), b, c, uind, uval)
Model(env, A::Ta, b, c, uind, uval) where{Tv<:Real, Ta<:Matrix{Tv}} = Model(env, sparse(float(A)), b, c, uind, uval)
Model() = Model(spzeros(0, 0), Vector{Float64}(0,), Vector{Float64}(0,))
Model(A, b, c) = Model(A, b, c, Vector{Int}(0,), Vector{Float64}(0,))


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
getnumvar(m::Model) = m.n_var

"""
    getnumconstr(m::Model)

Return the number of constraints in the model. This number does not include
    explicit bounds on the variables.
"""
getnumconstr(m::Model) = m.n_con

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
    size(c, 1) == m.n_var || throw(DimensionMismatch(
        "c has $(size(c, 1)) coeffs but model has $(m.n_var) variables"
    ))

    m.c = copy(c)
    return nothing
end

"""
    getvarlowerbounds(m::Model)

Return lower bounds on the variables.
"""
getvarlowerbounds(m::Model) = zeros(m.n_var)

"""
    getvarupperbounds(m::Model)

Return upper bounds on the variables. If a given variable has no explicit upper
    bound, the returned value is `Inf`.
"""
function getvarupperbounds(m::Model)
    ub = fill(Inf, m.n_var)
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
    size(ub, 1) == m.n_var || throw(DimensionMismatch(
        "ub has size $(size(c, 1)) but model has $(m.n_var) variables"
    ))
    minimum(ub) >= zero(Tv) || error("Upper bounds must be non-negative")

    m.uind = ub.nzind
    m.uval = ub.nzval
    m.n_var_ub = nnz(ub)
    return nothing
end
function setvarupperbounds!(m::Model, ub::AbstractArray{Tv}) where{Tv<:Real}
    # Dimension check
    size(ub, 1) == m.n_var || throw(DimensionMismatch(
        "ub has size $(size(c, 1)) but model has $(m.n_var) variables"
    ))
    minimum(ub) >= zero(Tv) || error("Upper bounds must be non-negative")

    u_ = ub .< Inf
    m.uind = collect(1:m.n_var)[u_]
    m.uval = ub[u_]
    m.n_var_ub = sum(u_)

    return nothing
end

"""
    addvar!(m::Model, colval, l, u, objcoeff)

Add a variable to the model.
"""
function addvar!(m::Model, colvals::AbstractVector{Tv}, l::Real, u::Real, objcoef::Real) where Tv<:Real

    # Dimension check
    m.n_con == size(colvals, 1) || throw(DimensionMismatch(
        "Column has $(size(col, 1)) coeffs but model has $(m.n_con) constraints"
    ))
    l <= u || error("Upper-bound must be greater than lower bound")

    if l == -Inf && u == Inf
        # free variable: add positive and negative parts
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, colvals)
        m.n_var += 1

        insert!(m.c, k, objcoef)

        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, -colvals)
        m.n_var += 1

        insert!(m.c, k, -objcoef)


    elseif l == -Inf && u < Inf
        # upperbounded, no lower bound
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, -colvals)
        m.n_var += 1

        # Update objective
        insert!(m.c, k, -objval)

        # Update right-hand side
        m.b .+= u * colvals

    elseif -Inf < l && u < Inf
        # lower- and upper-bounded
        m.A, k = Tulip.LinearAlgebra.addcolumn!(m.A, colvals)
        m.n_var += 1

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
        m.n_var += 1

        # Update objective
        insert!(m.c, k, objcoef)

        # Update right-hand side
        m.b .-= l * colvals
    end

    m.status = :Built

    return nothing
end
addvar!(m::Model, constridx, constrcoef, l, u, objcoef) = addvar!(m, sparsevec(constridx, constrcoef, m.n_con), l, u, objcoef)
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
function addconstr!(m::Model, rowvals::AbstractVector{Tv}, rhs::Real) where Tv<:Real

    # Dimension checks
    m.n_var == size(rowvals, 1) || throw(DimensionMismatch(
        "Row has $(size(rowals, 1)) coefs but model has $(m.n_var) variables"
    ))
    -Inf < rhs < Inf || error("Right-hand side must have finite value")

    # Add constraint
    m.n_con += 1
    m.b = vcat(m.b, rhs)
    m.A = vcat(m.A, rowvals)

    return nothing
end

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
getobjectivevalue(m::Model) = dot(m.c, m.x)

"""
    getdualbound(m::Model)

Return dual bound on the obective value. Returns a lower (resp. upper) bound
if the problem is a minimization (resp. maximization).
"""
getdualbound(m::Model) = dot(m.b, m.y) - dot(m.uval, m.z)

"""
    getobjectivedualgap(m::Model)

Return the duality gap. 
"""
getobjectivedualgap(m::Model) = dot(m.x, m.s) + dot(m.w, m.z)

"""
    getsolution(m::Model)

Return best known (primal) solution to the problem.
"""
getsolution(m::Model) = copy(m.x)

"""
    getconstrduals(m::Model)

Return dual variables associated to linear constraints.
"""
getconstrduals(m::Model) = copy(m.y)

"""
    getreducedcosts(m::Model)

Return reduced costs of primal variables.
"""
getreducedcosts(m::Model) = copy(m.s)

"""
    getnumbarrieriter(m::Model)

Return number of barrier iterations.
"""
getnumbarrieriter(m::Model) = m.numbarrieriter

"""
    getsolutiontime(m::Model)

Return runtime of the optimizer.
"""
getsolutiontime(m::Model) = m.runtime