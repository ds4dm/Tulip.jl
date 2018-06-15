import MathProgBase
const MPB = MathProgBase


export TulipSolver

"""
    TulipSolver

"""
struct TulipSolver <: MPB.AbstractMathProgSolver
    options
end

function TulipSolver(;kwargs...)

    TulipSolver(kwargs)
end

"""
    TulipMathProgModel

Wrapper for MathProgBase models.
"""
mutable struct TulipMathProgModel <: MPB.AbstractLinearQuadraticModel
    inner::Model
end


#=======================================================
    MathProgBase - SolverInterface
=======================================================#

function MPB.LinearQuadraticModel(s::TulipSolver)
    
    # create an empty model
    m = Tulip.Model()

    # set user-defined parameters
    for (name, val) in s.options
        m.env[name] = val
    end

    TulipMathProgModel(m)
end

MPB.getsolution(m::TulipMathProgModel) = copy(m.inner.sol.x)

MPB.getobjval(m::TulipMathProgModel) = dot(m.inner.c, m.inner.sol.x)

MPB.optimize!(m::TulipMathProgModel) = optimize!(m.inner)

MPB.status(m::TulipMathProgModel) = m.inner.status

function MPB.getobjbound(m::TulipMathProgModel)
    warn("Result may be wrong if current solution is not feasible.")
    return dot(m.inner.b, m.inner.sol.y) - dot(m.inner.uind, m.inner.sol.z)
end

function MPB.getobjgap(m::TulipMathProgModel)
    warn("Result may be wrong if current solution is not feasible.")
    return dot(m.inner.sol.x, m.inner.sol.s) + dot(m.inner.sol.w, m.inner.sol.z)
end

MPB.getrawsolver(m::TulipMathProgModel) = m.inner

function MPB.getsolvetime(m::TulipMathProgModel)
    warn("MPB.getsolvetime currently not implemented. Function call ignored.")
    return -1.0
end

function MPB.getsense(m::TulipMathProgModel)
    return :Min
end

function MPB.setsense!(m::TulipMathProgModel, sense)
    warn("MPB.setsense! currently not implemented. Function call ignored.")
    return nothing
end

MPB.numvar(m::TulipMathProgModel) = copy(m.inner.n_var)

MPB.numconstr(m::TulipMathProgModel) = copy(m.inner.n_con)




#=======================================================
    MathProgBase - LinearQuadratic
=======================================================#

function MPB.loadproblem!(
    m::TulipMathProgModel,
    A, l, u, c, lb, ub, sense
)
    
    n_var = size(A, 2)
    n_con = size(A, 1)

    # Dimension check
    n_var == size(c, 1) || throw(DimensionMismatch("c has size $(size(c))"))
    n_var == size(l, 1) || throw(DimensionMismatch("l has size $(size(c))"))
    n_var == size(u, 1) || throw(DimensionMismatch("u has size $(size(c))"))
    n_con == size(lb, 1) || throw(DimensionMismatch("lb has size $(size(c))"))
    n_con == size(ub, 1) || throw(DimensionMismatch("ub has size $(size(c))"))

    # Verify if problem is in standard form
    # TODO: handle conversion to standard form in model constructor
    if lb != ub
        error("Only equality constraints are supported.")
    end
    if sense != :Min
        error("Only minimization is supported.")
    end

    # extract upper bounds
    uind = Vector{Int}(0,)
    uval = Vector{Float64}(0,)
    n_var_ub = 0
    for i in 1:n_var
        if u[i] < Inf
            n_var_ub += 1
            push!(uind, i)
            push!(uval, u[i])
        end
    end

    m.inner = Model(copy(m.inner.env), A, ub, c, uind, uval)

    return nothing
end

MPB.getvarLB(m::TulipMathProgModel) = getvarlowerbounds(m.inner)

function MPB.setvarLB!(m::TulipMathProgModel, collb)
    warn("MPB.setvarLB! currently not implemented. Function call ignored.")
    return nothing
end

MPB.getvarUB(m::TulipMathProgModel) = getvarupperbounds(m.inner)

MPB.setvarUB!(m::TulipMathProgModel, colub) = setvarupperbounds!(m.inner, colub)

MPB.getconstrLB(m::TulipMathProgModel) = getconstrlowerbound(m.inner)

function MPB.setconstrLB!(m::TulipMathProgModel, rowlb)
    warn("MPB.setconstrLB! currently not implemented. Function call ignored.")
    return nothing
end

MPB.getconstrUB(m::TulipMathProgModel) = getconstrupperbound(m.inner)

function MPB.setconstrUB!(m::TulipMathProgModel, rowub)
    warn("MPB.setconstrUB! currently not implemented. Function call ignored.")
    return nothing
end

MPB.getobj(m::TulipMathProgModel) = getobjectivecoeffs(m.inner)

MPB.setobj!(m::TulipMathProgModel, c) = setobjectivecoeffs!(m, c)

MPB.getconstrmatrix(m::TulipMathProgModel) = getlinearconstrcoeffs(m.inner)

MPB.addvar!(m::TulipMathProgModel, constridx, constrcoef, l, u, objcoef) = addvar!(m.inner, constridx, constrcoef, l, u, objcoef)
MPB.addvar!(m::TulipMathProgModel, l, u, objcoef) = MPB.addvar!(m, [], [], l, u, objcoef)

function MPB.delvars!(m::TulipMathProgModel, idxs)
    warn("MPB.delvars! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.addconstr!(m::TulipMathProgModel, varidx, coef, lb, ub)
    warn("MPB.addconstrs! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.delconstrs!(m::TulipMathProgModel, idxs)
    warn("MPB.delconstrs! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.changecoeffs!(m::TulipMathProgModel, cidxs, vidxs, val)
    warn("MPB.changecoeffs! currently not implemented. Function call ignored.")
    return nothing
end

MPB.numlinconstr(m::TulipMathProgModel) = getnumconstr(m.inner)

function MPB.getconstrsolution(m::TulipMathProgModel)
    return copy(m.inner.A * m.inner.sol.x)
end

MPB.getreducedcosts(m::TulipMathProgModel) = getreducedcosts(m.inner)

MPB.getconstrduals(m::TulipMathProgModel) = getconstrduals(m.inner)

function MPB.getinfeasibilityray(m::TulipMathProgModel)
    warn("MPB.getinfeasibilityray! currently not implemented. Function call ignored.")
    return Vector{Float64}(0,)
end

function MPB.getbarrieriter(m::TulipMathProgModel)
    warn("MPB.getbarrieriter! currently not implemented. Function call ignored.")
    return -1
end