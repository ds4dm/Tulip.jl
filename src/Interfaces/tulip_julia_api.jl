using QPSReader

# TODO: user-facing API in Julia
# Other APIs should wrap this one
# TODO: docstrings
# TODO: define Traits on attributes (e.g.: IsModifiable, IsNumeric, etc..)
#   for error messages
abstract type AbstractAttribute end

# Model attributes
struct ModelName          <: AbstractAttribute end
struct NumberOfVariables  <: AbstractAttribute end
struct ObjectiveValue     <: AbstractAttribute end
struct DualObjectiveValue <: AbstractAttribute end
struct ObjectiveConstant  <: AbstractAttribute end
struct ObjectiveSense     <: AbstractAttribute end
struct Status             <: AbstractAttribute end
struct BarrierIterations  <: AbstractAttribute end
struct SolutionTime       <: AbstractAttribute end

# Variable attributes
struct VariableLowerBound <: AbstractAttribute end
struct VariableUpperBound <: AbstractAttribute end
struct VariableObjective  <: AbstractAttribute end
struct VariableName       <: AbstractAttribute end

# Constraint attributes
struct ConstraintLowerBound <: AbstractAttribute end
struct ConstraintUpperBound <: AbstractAttribute end
struct ConstraintName       <: AbstractAttribute end

# Model creation and modification
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

function set_attribute(m::Model{Tv}, ::VariableLowerBound, j::Int, lb::Tv) where{Tv}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.lvar[j] = lb
    return nothing
end

function set_attribute(m::Model{Tv}, ::VariableUpperBound, j::Int, ub::Tv) where{Tv}
    # sanity checks
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")

    # Update bound
    m.pbdata.uvar[j] = ub
    return nothing
end

function set_attribute(m::Model{Tv}, ::ConstraintLowerBound, i::Int, lb::Tv) where{Tv}
    # sanity checks
    1 <= i <= m.pbdata.nvar || error("Invalid constraint index $i")

    # Update bound
    m.pbdata.lcon[i] = lb
    return nothing
end

function set_attribute(m::Model{Tv}, ::ConstraintUpperBound, i::Int, ub::Tv) where{Tv}
    # sanity checks
    1 <= i <= m.pbdata.nvar || error("Invalid constraint index $i")

    # Update bound
    m.pbdata.ucon[i] = ub
    return nothing
end

function get_attribute(m::Model, ::VariableName, j::Int)
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")
    return m.pbdata.var_names[j]
end

function set_attribute(m::Model, ::VariableName, j::Int, name::String)
    1 <= j <= m.pbdata.nvar || error("Invalid variable index $j")
    # TODO: ensure that no two variables have the same name
    m.pbdata.var_names[j] = name
    return nothing
end

function get_attribute(m::Model, ::ConstraintName, i::Int)
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")
    return m.pbdata.con_names[i]
end

function set_attribute(m::Model, ::ConstraintName, i::Int, name::String)
    1 <= i <= m.pbdata.ncon || error("Invalid constraint index $i")
    m.pbdata.con_names[i] = name
    return nothing
end

# TODO: Set/get parameters

# TODO: Query solution value and attributes
get_attribute(m::Model, ::ObjectiveConstant) = m.pbdata.obj0
set_attribute(m::Model{Tv}, ::ObjectiveConstant, obj0::Tv) where{Tv} = (m.pbdata.obj0 = obj0; return nothing)


function get_attribute(m::Model{Tv}, ::ObjectiveValue) where{Tv}
    if isnothing(m.solver)
        error("No solver is attached to the model")
    else
        z = m.solver.primal_bound_scaled
        return m.pbdata.objsense ? z : -z
    end
end

function get_attribute(m::Model{Tv}, ::DualObjectiveValue) where{Tv}
    if isnothing(m.solver)
        error("No solver is attached to the model")
    else
        z = m.solver.dual_bound_scaled
        return m.pbdata.objsense ? z : -z
    end
end