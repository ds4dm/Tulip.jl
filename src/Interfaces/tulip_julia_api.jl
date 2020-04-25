using QPSReader

# TODO: user-facing API in Julia
# Other APIs should wrap this one
# TODO: docstrings
# TODO: define Traits on attributes (e.g.: IsModifiable, IsNumeric, etc..)
#   for error messages


"""
    load_problem!(m::Model{Tv}, fname::String)

Read a model from file `fname` and load it into model `m`.

Only free MPS files are currently supported.
"""
function load_problem!(m::Model{Tv}, fname::String) where{Tv}
    Base.empty!(m)

    dat = with_logger(Logging.NullLogger()) do
        readqps(fname, mpsformat=:free)
    end

    # TODO: avoid allocations when Tv is Float64
    objsense = !(dat.objsense == :max)
    load_problem!(m.pbdata,
        dat.name,
        objsense, Tv.(dat.c), Tv(dat.c0),
        sparse(dat.arows, dat.acols, Tv.(dat.avals), dat.ncon, dat.nvar),
        Tv.(dat.lcon), Tv.(dat.ucon),
        Tv.(dat.lvar), Tv.(dat.uvar),
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
    set_attribute(m::Model{Tv}, ::VariableLowerBound, j::Int, lb::Tv)

Set the lower bound of variable `j` in model `m` to `lb`.
"""
function set_attribute(m::Model{Tv}, ::VariableLowerBound, j::Int, lb::Tv) where{Tv}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.lvar[j] = lb
    return nothing
end

"""
    get_attribute(m::Model{Tv}, ::VariableLowerBound, j::Int)

Query the lower bound of variable `j` in model `m`.
"""
function get_attribute(m::Model, ::VariableLowerBound, j::Int)
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    return m.pbdata.lvar[j]
end

"""
    set_attribute(m::Model{Tv}, ::VariableUpperBound, j::Int, ub::Tv)

Set the upper bound of variable `j` in model `m` to `ub`.
"""
function set_attribute(m::Model{Tv}, ::VariableUpperBound, j::Int, ub::Tv) where{Tv}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.uvar[j] = ub
    return nothing
end

"""
    set_attribute(m::Model{Tv}, ::ConstraintLowerBound, i::Int, lb::Tv)

Set the lower bound of constraint `i` in model `m` to `lb`.
"""
function set_attribute(m::Model{Tv}, ::ConstraintLowerBound, i::Int, lb::Tv) where{Tv}
    # sanity checks
    1 <= i <= m.pbdata.nvar || error("Invalid constraint index $i")

    # Update bound
    m.pbdata.lcon[i] = lb
    return nothing
end

"""
    set_attribute(m::Model{Tv}, ::ConstraintUpperBound, i::Int, ub::Tv)

Set the upper bound of constraint `i` in model `m` to `ub`.
"""
function set_attribute(m::Model{Tv}, ::ConstraintUpperBound, i::Int, ub::Tv) where{Tv}
    # sanity checks
    1 <= i <= m.pbdata.nvar || error("Invalid constraint index $i")

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
    setfield!(m.params, Symbol(pname), val)
    return nothing
end


# TODO: Query solution value
get_attribute(m::Model, ::ObjectiveConstant) = m.pbdata.obj0
set_attribute(m::Model{Tv}, ::ObjectiveConstant, obj0::Tv) where{Tv} = (m.pbdata.obj0 = obj0; return nothing)

"""
    get_attribute(model::Model, ::ObjectiveValue)

Query the `ObjectiveValue` attribute from `model`
"""
function get_attribute(m::Model{Tv}, ::ObjectiveValue) where{Tv}
    if isnothing(m.solution)
        error("Model has no solution")
    else
        z = m.solution.z_primal
        return m.pbdata.objsense ? z : -z
    end
end

"""
    get_attribute(model::Model, ::DualObjectiveValue)

Query the `DualObjectiveValue` attribute from `model`
"""
function get_attribute(m::Model{Tv}, ::DualObjectiveValue) where{Tv}
    if isnothing(m.solution)
        error("Model has no solution")
    else
        z = m.solution.z_dual
        return m.pbdata.objsense ? z : -z
    end
end