import MathProgBase
const MPB = MathProgBase


export TulipSolver

"""
    TulipSolver

"""
struct TulipSolver <: MPB.AbstractMathProgSolver
    #=======================================================
        Optimization-related parameters
    =======================================================#

    # I/O behaviour
    output_level::Int       # 0 means no output, 1 means normal

    # Termination criteria
    n_iter_max::Int         # Maximum number of barrier iterations
    time_limit::Float64     # Time limit (in seconds)
    eps_tol_p::Float64      # Numerical tolerance for primal feasibility
    eps_tol_d::Float64      # Numerical tolerance for dual feasibility
    eps_tol_g::Float64      # Numerical tolerance for optimality gap
end

function TulipSolver(;
    output_level=1,
    n_iter_max=100,
    time_limit=Inf,
    eps_tol_p=10.0^-8,
    eps_tol_d=10.0^-8,
    eps_tol_g=10.0^-8
)

    TulipSolver(output_level, n_iter_max, time_limit, eps_tol_p, eps_tol_d, eps_tol_g)
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

    inner_model = Tulip.Model(
        output_level=s.output_level,
        n_iter_max=s.n_iter_max,
        time_limit=s.time_limit,
        eps_tol_p=s.eps_tol_p,
        eps_tol_d=s.eps_tol_d,
        eps_tol_g=s.eps_tol_g
    )
    
    TulipMathProgModel(inner_model)
end

function MPB.getsolution(m::TulipMathProgModel)
    return copy(m.inner.sol.x)
end

function MPB.getobjval(m::TulipMathProgModel)
    return dot(m.inner.c, m.inner.sol.x)
end

function MPB.optimize!(m::TulipMathProgModel)

    optimize!(m.inner)

    return m.inner.status

end

function MPB.status(m::TulipMathProgModel)
    return copy(m.inner.status)
end

function MPB.getobjbound(m::TulipMathProgModel)
    warn("Result may be wrong if current solution is not feasible.")
    return dot(m.inner.b, m.inner.sol.y) - dot(m.inner.uind, m.inner.sol.z)
end

function MPB.getobjgap(m::Tulip)
    warn("Result may be wrong if current solution is not feasible.")
    return dot(m.inner.sol.x, m.inner.sol.s) + dot(m.inner.sol.w, m.inner.sol.z)
end

function MPB.getrawsolver(m::TulipMathProgModel)
    return m.inner
end

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

function MPB.numvar(m::TulipMathProgModel)
    return copy(m.inner.n_var)
end

function MPB.numconstr(m::TulipMathProgModel)
    return copy(m.inner.n_con)
end




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

    m.inner.n_var = n_var
    m.inner.n_con = n_con
    m.inner.n_var_ub = n_var_ub
    m.inner.A = A
    m.inner.b = ub
    m.inner.c = c
    m.inner.uind = uind
    m.inner.uval = uval

    m.inner.sol = PrimalDualPoint(
        ones(n_var),
        ones(n_var_ub),
        zeros(n_con),
        ones(n_var),
        ones(n_var_ub)
    )

    m.inner.status = :Built

    return nothing
end

function MPB.getvarLB(m::TulipMathProgModel)
    return zeros(m.inner.n_var)
end

function MPB.setvarLB!(m::TulipMathProgModel, collb)
    warn("MPB.setvarLB! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.getvarUB(m::TulipMathProgModel)

    u = Inf*ones(m.inner.n_var)
    u[m.inner.uind] = m.inner.uval
    
    return u
end

function MPB.setvarUB!(m::TulipMathProgModel, colub)
    warn("MPB.setvarUB! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.getconstrLB(m::TulipMathProgModel)
    return copy(m.inner.b)
end

function MPB.setconstrLB!(m::TulipMathProgModel, rowlb)
    warn("MPB.setconstrLB! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.getconstrUB(m::TulipMathProgModel)
    return copy(m.inner.b)
end

function MPB.setconstrUB!(m::TulipMathProgModel, rowub)
    warn("MPB.setconstrUB! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.getobj(m::TulipMathProgModel)
    return copy(m.inner.c)
end

function MPB.setobj!(m::TulipMathProgModel, c)
    warn("MPB.setobj! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.getconstrmatrix(m::TulipMathProgModel)
    return copy(m.inner.A)
end

function MPB.addvar!(m::TulipMathProgModel, l, u, objcoef)
    MPB.addvar!(m, [], [], l, u, objcoef)
end

function MPB.addvar!(m::TulipMathProgModel, constridx, constrcoef, l, u, objcoef)
    warn("MPB.addvar! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.delvars!(m::TulipMathProgModel, idxs)
    warn("MPB.delvars! currently not implemented. Function call ignored.")
    return nothing
end

function MPB.addconstrs!(m::TulipMathProgModel, varidx, coef, lb, ub)
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

function MPB.numlinconstr(m::TulipMathProgModel)
    return copy(m.inner.n_con)
end

function MPB.getconstrsolution(m::TulipMathProgModel)
    return copy(m.inner.A * m.inner.sol.x)
end

function MPB.getreducedcosts(m::TulipMathProgModel)
    warn("MPB.changecoeffs! currently not implemented. Function call ignored.")
    return Vector{Float64}(0,)
end

function MPB.getconstrduals(m::TulipMathProgModel)
    return copy(m.inner.sol.y)
end

function MPB.getinfeasibilityray(m::TulipMathProgModel)
    warn("MPB.getinfeasibilityray! currently not implemented. Function call ignored.")
    return Vector{Float64}(0,)
end

function MPB.getbarrieriter(m::TulipMathProgModel)
    warn("MPB.getbarrieriter! currently not implemented. Function call ignored.")
    return -1
end