# =============================================
#   Supported attributes
# =============================================
const SUPPORTED_OPTIMIZER_ATTR = Union{
    MOI.NumberOfThreads,
    MOI.RawOptimizerAttribute,
    MOI.SolverName,
    MOI.SolverVersion,
    MOI.SolveTimeSec,
    MOI.Silent,
    MOI.TimeLimitSec,
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_OPTIMIZER_ATTR} = true


# =============================================
#   1. Optimizer attributes
# =============================================

#
#   NumberOfThreads
#
MOI.get(m::Optimizer, ::MOI.NumberOfThreads) = m.inner.params.Threads

function MOI.set(m::Optimizer, ::MOI.NumberOfThreads, n::Int)
    # TODO: use lower-level API
    m.inner.params.Threads = n
    return nothing
end

#
#   SolverName
#
MOI.get(::Optimizer, ::MOI.SolverName) = "Tulip"

#
#   SolverVersion
#
MOI.get(::Optimizer, ::MOI.SolverVersion) = string(Tulip.version())

#
#   SolveTimeSec
#
MOI.get(m::Optimizer, ::MOI.SolveTimeSec) = m.solve_time

#
#   Silent
#
MOI.get(m::Optimizer, ::MOI.Silent) = m.inner.params.OutputLevel <= 0

function MOI.set(m::Optimizer, ::MOI.Silent, flag::Bool)
    m.inner.params.OutputLevel = 1 - flag
    # TODO: make a decision about LogLevel
    return nothing
end

#
#   TimeLimitSec
#
function MOI.get(m::Optimizer, ::MOI.TimeLimitSec)
    value = m.inner.params.IPM.TimeLimit
    return value == Inf ? nothing : value
end

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, t::Union{Real,Nothing})
    m.inner.params.IPM.TimeLimit = convert(Float64, something(t, Inf))
    return nothing
end

#
#   RawParameter
#
MOI.get(m::Optimizer, attr::MOI.RawOptimizerAttribute) = get_parameter(m.inner, attr.name)

MOI.set(m::Optimizer, attr::MOI.RawOptimizerAttribute, val) = set_parameter(m.inner, attr.name, val)


# =============================================
#   2. Model attributes
# =============================================
const SUPPORTED_MODEL_ATTR = Union{
    MOI.Name,
    MOI.ObjectiveSense,
    MOI.NumberOfVariables,
    MOI.ListOfVariableIndices,
    MOI.ListOfConstraintIndices,
    MOI.NumberOfConstraints,
    # ListOfConstraints,  # TODO
    MOI.ObjectiveFunctionType,
    MOI.ObjectiveValue,
    MOI.DualObjectiveValue,
    MOI.RelativeGap,
    # MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    MOI.RawStatusString,
    MOI.ResultCount,
    MOI.TerminationStatus,
    MOI.PrimalStatus,
    MOI.DualStatus
}

MOI.supports(::Optimizer, ::SUPPORTED_MODEL_ATTR) = true

#
#   ListOfModelAttributesSet
#
function MOI.get(m::Optimizer{T}, ::MOI.ListOfModelAttributesSet) where {T}
    ret = MOI.AbstractModelAttribute[]
    if !isempty(m.inner.pbdata.name)
        push!(ret, MOI.Name())
    end
    if m.objective_sense !== nothing
        push!(ret, MOI.ObjectiveSense())
    end
    if m._obj_type == _SINGLE_VARIABLE
        push!(ret, MOI.ObjectiveFunction{MOI.VariableIndex}())
    elseif m._obj_type == _SCALAR_AFFINE
        push!(ret, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
    end
    return ret
end

#
#   ListOfVariableIndices
#
function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return copy(m.var_indices_moi)
end

#
#   Name
#
MOI.get(m::Optimizer, ::MOI.Name) = m.inner.pbdata.name
MOI.set(m::Optimizer, ::MOI.Name, name) = (m.inner.pbdata.name = name)

#
#   NumberOfVariables
#
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = m.inner.pbdata.nvar

#
#   ObjectiveFunctionType
#
function MOI.get(m::Optimizer{T}, ::MOI.ObjectiveFunctionType) where{T}
    if m._obj_type == _SINGLE_VARIABLE
        return MOI.VariableIndex
    end
    return MOI.ScalarAffineFunction{T}
end

#
#   ObjectiveSense
#
function MOI.get(m::Optimizer, ::MOI.ObjectiveSense)
    return something(m.objective_sense, MOI.FEASIBILITY_SENSE)
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense)
    m.objective_sense = s
    if s == MOI.MIN_SENSE || s == MOI.FEASIBILITY_SENSE
        m.inner.pbdata.objsense = true
    else
        @assert s == MOI.MAX_SENSE
        m.inner.pbdata.objsense = false
    end
    return nothing
end

#
#   ObjectiveValue
#
function MOI.get(m::Optimizer{T}, attr::MOI.ObjectiveValue) where{T}
    MOI.check_result_index_bounds(m, attr)
    raw_z = get_attribute(m.inner, ObjectiveValue())
    is_feas = MOI.get(m, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    return raw_z * !is_feas
end

#
#   DualObjectiveValue
#
function MOI.get(m::Optimizer{T}, attr::MOI.DualObjectiveValue) where{T}
    MOI.check_result_index_bounds(m, attr)
    return get_attribute(m.inner, DualObjectiveValue())
end

MOI.get(m::Optimizer, ::MOI.ObjectiveBound) = MOI.get(m, MOI.DualObjectiveValue())

#
#   RawSolver
#
MOI.get(m::Optimizer, ::MOI.RawSolver) = m.inner

#
#   RelativeGap
#
function MOI.get(m::Optimizer{T}, ::MOI.RelativeGap) where{T}
    # TODO: dispatch a function call on m.inner
    zp = m.inner.solver.primal_objective
    zd = m.inner.solver.dual_objective
    return (abs(zp - zd) / (T(1 // 10^6)) + abs(zd))
end

#
#   RawStatusString
#
function MOI.get(m::Optimizer, ::MOI.RawStatusString)
    return string(m.inner.status)
end

#
#   ResultCount
#
function MOI.get(m::Optimizer, ::MOI.ResultCount)
    st = MOI.get(m, MOI.TerminationStatus())

    if (st == MOI.OPTIMIZE_NOT_CALLED
        || st == MOI.OTHER_ERROR
        || st == MOI.MEMORY_LIMIT
    )
        return 0
    end
    return 1
end

#
#   SimplexIterations
#
MOI.get(::Optimizer, ::MOI.SimplexIterations) = 0

#
#   BarrierIterations
#
# TODO: use inner query
MOI.get(m::Optimizer, ::MOI.BarrierIterations) = get_attribute(m.inner, BarrierIterations())

#
#   TerminationStatus
#
# TODO: use inner query
function MOI.get(m::Optimizer, ::MOI.TerminationStatus)
    return MOITerminationStatus(get_attribute(m.inner, Status()))
end

#
#   PrimalStatus
#
# TODO: use inner query
function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    attr.result_index == 1 || return MOI.NO_SOLUTION

    if isnothing(m.inner.solution)
        return MOI.NO_SOLUTION
    else
        MOISolutionStatus(m.inner.solution.primal_status)
    end
end

#
#   DualStatus
#
# TODO: use inner query
function MOI.get(m::Optimizer, attr::MOI.DualStatus)
    attr.result_index == 1 || return MOI.NO_SOLUTION

    if isnothing(m.inner.solution)
        return MOI.NO_SOLUTION
    else
        MOISolutionStatus(m.inner.solution.dual_status)
    end
end
