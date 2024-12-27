using QPSReader
using TimerOutputs: tottime

# TODO: user-facing API in Julia
# Other APIs should wrap this one
# TODO: docstrings
# TODO: define Traits on attributes (e.g.: IsModifiable, IsNumeric, etc..)
#   for error messages


"""
    load_problem!(m::Model{T}, fname::String)

Read a model from file `fname` and load it into model `m`.

Only free MPS files are currently supported.
"""
function load_problem!(m::Model{T}, fname::String) where{T}
    Base.empty!(m)

    dat = with_logger(Logging.NullLogger()) do
        _open(fname) do io
            readqps(io, mpsformat=:free)
        end
    end

    # TODO: avoid allocations when T is Float64
    objsense = !(dat.objsense == :max)
    load_problem!(m.pbdata,
        dat.name,
        objsense, T.(dat.c), T(dat.c0),
        sparse(dat.arows, dat.acols, T.(dat.avals), dat.ncon, dat.nvar),
        T.(dat.lcon), T.(dat.ucon),
        T.(dat.lvar), T.(dat.uvar),
        dat.connames, dat.varnames
    )

    return m
end

"""
    get_attribute(model::Model, ::ModelName)

Query the `ModelName` attribute from `model`
"""
get_attribute(m::Model, ::ModelName) = m.pbdata.name

"""
    set_attribute(model::Model, ::ModelName, name::String)

Set the `ModelName` attribute in `model`
"""
set_attribute(m::Model, ::ModelName, name::String) = (m.pbdata.name = name; return nothing)

"""
    get_attribute(model::Model, ::Status)

Query the `Status` attribute from `model`
"""
function get_attribute(m::Model, ::Status)
    return m.status
end

"""
    get_attribute(model::Model, ::SolutionTime)

Query the `SolutionTime` attribute from `model`
"""
function get_attribute(m::Model, ::SolutionTime)
    if isnothing(m.solver)
        return 0
    else
      local ns = tottime(m.solver.timer)
      return ns * 1e-9
    end
end

"""
    get_attribute(model::Model, ::BarrierIterations)

Query the `BarrierIterations` attribute from `model`
"""
function get_attribute(m::Model, ::BarrierIterations)
    if isnothing(m.solver)
        return 0
    else
        return m.solver.niter
    end
end

"""
    set_attribute(m::Model{T}, ::VariableLowerBound, j::Int, lb::T)

Set the lower bound of variable `j` in model `m` to `lb`.
"""
function set_attribute(m::Model{T}, ::VariableLowerBound, j::Int, lb::T) where{T}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.lvar[j] = lb
    return nothing
end

"""
    get_attribute(m::Model{T}, ::VariableLowerBound, j::Int)

Query the lower bound of variable `j` in model `m`.
"""
function get_attribute(m::Model, ::VariableLowerBound, j::Int)
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    return m.pbdata.lvar[j]
end

"""
    set_attribute(m::Model{T}, ::VariableUpperBound, j::Int, ub::T)

Set the upper bound of variable `j` in model `m` to `ub`.
"""
function set_attribute(m::Model{T}, ::VariableUpperBound, j::Int, ub::T) where{T}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.uvar[j] = ub
    return nothing
end

"""
    set_attribute(m::Model{T}, ::ConstraintLowerBound, i::Int, lb::T)

Set the lower bound of constraint `i` in model `m` to `lb`.
"""
function set_attribute(m::Model{T}, ::ConstraintLowerBound, i::Int, lb::T) where{T}
    # sanity checks
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")

    # Update bound
    m.pbdata.lcon[i] = lb
    return nothing
end

"""
    set_attribute(m::Model{T}, ::ConstraintUpperBound, i::Int, ub::T)

Set the upper bound of constraint `i` in model `m` to `ub`.
"""
function set_attribute(m::Model{T}, ::ConstraintUpperBound, i::Int, ub::T) where{T}
    # sanity checks
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")

    # Update bound
    m.pbdata.ucon[i] = ub
    return nothing
end

"""
    get_attribute(m::Model, ::VariableName, j::Int)

Query the name of variable `j` in model `m`
"""
function get_attribute(m::Model, ::VariableName, j::Int)
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")
    return m.pbdata.var_names[j]
end

"""
    set_attribute(m::Model, ::VariableName, j::Int, name::String)

Set the name of variable `j` in model `m` to `name`.
"""
function set_attribute(m::Model, ::VariableName, j::Int, name::String)
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")
    # TODO: ensure that no two variables have the same name
    m.pbdata.var_names[j] = name
    return nothing
end

"""
    get_attribute(m::Model, ::ConstraintName, i::Int)

Query the name of constraint `i` in model `m`
"""
function get_attribute(m::Model, ::ConstraintName, i::Int)
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")
    return m.pbdata.con_names[i]
end

"""
    set_attribute(m::Model, ::ConstraintName, i::Int, name::String)

Set the name of constraint `i` in model `m` to `name`.
"""
function set_attribute(m::Model, ::ConstraintName, i::Int, name::String)
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")
    m.pbdata.con_names[i] = name
    return nothing
end

# TODO: Set/get parameters
"""
    get_parameter(m::Model, pname::String)

Query the value of parameter `pname` in model `m`.
"""
function get_parameter(m::Model, pname::String)
    return getfield(m.params, Symbol(pname))
end

"""
    set_parameter(m::Model, pname::String, val)

Set the value of parameter `pname` in model `m` to `val`
"""
function set_parameter(m::Model, pname::String, val)
    if length(pname) > 4 && pname[1:4] == "IPM_"
        setfield!(m.params.IPM, Symbol(pname[5:end]), val)
    elseif length(pname) > 4 && pname[1:4] == "KKT_"
        setfield!(m.params.KKT, Symbol(pname[5:end]), val)
    elseif length(pname) > 10 && pname[1:9] == "Presolve_"
        setfield!(m.params.Presolve, Symbol(pname[10:end]), val)
    elseif hasfield(typeof(m.params), Symbol(pname))
        setfield!(m.params, Symbol(pname), val)
    else
        error("Unknown option: $pname")
    end
    return nothing
end


# TODO: Query solution value
get_attribute(m::Model, ::ObjectiveConstant) = m.pbdata.obj0
set_attribute(m::Model{T}, ::ObjectiveConstant, obj0::T) where{T} = (m.pbdata.obj0 = obj0; return nothing)

"""
    get_attribute(model::Model, ::ObjectiveValue)

Query the `ObjectiveValue` attribute from `model`
"""
function get_attribute(m::Model{T}, ::ObjectiveValue) where{T}
    if isnothing(m.solution)
        error("Model has no solution")
    end

    pst = m.solution.primal_status

    if pst != Sln_Unknown
        z  = dot(m.solution.x, m.pbdata.obj)
        # If solution is a ray, ignore constant objective term
        is_ray = m.solution.is_primal_ray
        z0 = !is_ray * m.pbdata.obj0
        return (z + z0)
    else
        # No solution, return zero
        return zero(T)
    end
end

"""
    get_attribute(model::Model, ::DualObjectiveValue)

Query the `DualObjectiveValue` attribute from `model`
"""
function get_attribute(m::Model{T}, ::DualObjectiveValue) where{T}
    if isnothing(m.solution)
        error("Model has no solution")
    end

    dst = m.solution.dual_status

    if dst != Sln_Unknown
        yl = m.solution.y_lower
        yu = m.solution.y_upper
        sl = m.solution.s_lower
        su = m.solution.s_upper

        bl = m.pbdata.lcon
        bu = m.pbdata.ucon
        xl = m.pbdata.lvar
        xu = m.pbdata.uvar

        z = (
              dot(yl, Diagonal(isfinite.(bl)), bl)
            - dot(yu, Diagonal(isfinite.(bu)), bu)
            + dot(sl, Diagonal(isfinite.(xl)), xl)
            - dot(su, Diagonal(isfinite.(xu)), xu)
        )

        # If solution is a ray, ignore constant objective term
        is_ray = m.solution.is_dual_ray
        z0 = !is_ray * m.pbdata.obj0
        return (z + z0)
    else
        # No solution, return zero
        return zero(T)
    end
end
