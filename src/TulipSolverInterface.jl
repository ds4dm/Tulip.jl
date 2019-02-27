import MathProgBase
const MPB = MathProgBase


export TulipSolver

"""
    TulipSolver

"""
struct TulipSolver <: MPB.AbstractMathProgSolver
    options
end

TulipSolver(;kwargs...) = TulipSolver(kwargs)

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
    env = TulipEnv()
    # set user-defined parameters
    for (name, val) in s.options
        env[name] = val
    end

    # Instanciate model
    m = Tulip.Model(env)
    TulipMathProgModel(m)
end

MPB.getsolution(m::TulipMathProgModel) = getprimalsolution(m.inner)

MPB.getobjval(m::TulipMathProgModel) = getobjectivevalue(m.inner)

MPB.optimize!(m::TulipMathProgModel) = optimize!(m.inner)

function MPB.status(m::TulipMathProgModel)
    
    s = m.inner.sln_status

    if s == Sln_Optimal
        return :Optimal
    elseif s == Sln_PrimalInfeasible || s == Sln_PrimalDualInfeasible
        return :Infeasible
    elseif s == Sln_DualInfeasible
        return :Unbounded
    else
        return :Unknown
    end

end

MPB.getobjbound(m::TulipMathProgModel) = getdualbound(m.inner)

MPB.getobjgap(m::TulipMathProgModel) = getobjectivedualgap(m.inner)

MPB.getrawsolver(m::TulipMathProgModel) = m.inner

MPB.getsolvetime(m::TulipMathProgModel) = m.inner.time_total

MPB.getsense(m::TulipMathProgModel) = :Min

function MPB.setsense!(m::TulipMathProgModel, sense)
    @warn "MPB.setsense! currently not implemented. Function call ignored."
    return nothing
end

MPB.numvar(m::TulipMathProgModel) = getnumvar(m.inner)

MPB.numconstr(m::TulipMathProgModel) = getnumconstr(m.inner)




#=======================================================
    MathProgBase - LinearQuadratic
=======================================================#

function MPB.loadproblem!(
    m::TulipMathProgModel,
    A, l, u, c, lb, ub, sense
)
    sense == :Min || error("Only minimization is supported.")
    
    loadmodel!(m.inner, size(A, 2), size(A, 1), A, c, l, u, lb, ub)
end

MPB.getvarLB(m::TulipMathProgModel) = getvarlowerbounds(m.inner)

function MPB.setvarLB!(m::TulipMathProgModel, collb)
    @warn "MPB.setvarLB! currently not implemented. Function call ignored."
    return nothing
end

MPB.getvarUB(m::TulipMathProgModel) = getvarupperbounds(m.inner)

MPB.setvarUB!(m::TulipMathProgModel, colub) = setvarupperbounds!(m.inner, colub)

MPB.getconstrLB(m::TulipMathProgModel) = getconstrlowerbounds(m.inner)

function MPB.setconstrLB!(m::TulipMathProgModel, rowlb)
    @warn "MPB.setconstrLB! currently not implemented. Function call ignored."
    return nothing
end

MPB.getconstrUB(m::TulipMathProgModel) = getconstrupperbounds(m.inner)

function MPB.setconstrUB!(m::TulipMathProgModel, rowub)
    @warn "MPB.setconstrUB! currently not implemented. Function call ignored."
    return nothing
end

MPB.getobj(m::TulipMathProgModel) = getobjectivecoeffs(m.inner)

MPB.setobj!(m::TulipMathProgModel, c) = setobjectivecoeffs!(m, c)

MPB.getconstrmatrix(m::TulipMathProgModel) = getlinearconstrcoeffs(m.inner)

MPB.addvar!(m::TulipMathProgModel, constridx, constrcoef, l, u, objcoef) = addvar!(m.inner, constridx, constrcoef, l, u, objcoef)
MPB.addvar!(m::TulipMathProgModel, l, u, objcoef) = MPB.addvar!(m, [], [], l, u, objcoef)

function MPB.delvars!(m::TulipMathProgModel, idxs)
    @warn "MPB.delvars! currently not implemented. Function call ignored."
    return nothing
end

function MPB.addconstr!(m::TulipMathProgModel, varidx, coef, lb, ub)
    @warn "MPB.addconstrs! currently not implemented. Function call ignored."
    return nothing
end

function MPB.delconstrs!(m::TulipMathProgModel, idxs)
    @warn "MPB.delconstrs! currently not implemented. Function call ignored."
    return nothing
end

function MPB.changecoeffs!(m::TulipMathProgModel, cidxs, vidxs, val)
    @warn "MPB.changecoeffs! currently not implemented. Function call ignored."
    return nothing
end

MPB.numlinconstr(m::TulipMathProgModel) = getnumconstr(m.inner)

MPB.getconstrsolution(m::TulipMathProgModel) = copy(m.inner.A * m.inner.sol.x)

MPB.getreducedcosts(m::TulipMathProgModel) = getreducedcosts(m.inner)

MPB.getconstrduals(m::TulipMathProgModel) = getdualsolution(m.inner)

MPB.getinfeasibilityray(m::TulipMathProgModel) = getinfeasibilityray(m.inner)

MPB.getunboundedray(m::TulipMathProgModel) = getunboundedray(m.inner)

MPB.getbarrieriter(m::TulipMathProgModel) = getnumbarrieriter(m.inner)