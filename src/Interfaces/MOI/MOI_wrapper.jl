import MathOptInterface as MOI

# ==============================================================================
#           HELPER FUNCTIONS
# ==============================================================================

"""
    MOITerminationStatus(st::TerminationStatus)

Convert a Tulip `TerminationStatus` into a `MOI.TerminationStatusCode`.
"""
function MOITerminationStatus(st::TerminationStatus)::MOI.TerminationStatusCode
    if st == Trm_NotCalled
        return MOI.OPTIMIZE_NOT_CALLED
    elseif st == Trm_Optimal
        return MOI.OPTIMAL
    elseif st == Trm_PrimalInfeasible
        return MOI.INFEASIBLE
    elseif st == Trm_DualInfeasible
        return MOI.DUAL_INFEASIBLE
    elseif st == Trm_IterationLimit
        return MOI.ITERATION_LIMIT
    elseif st == Trm_TimeLimit
        return MOI.TIME_LIMIT
    elseif st == Trm_MemoryLimit
        return MOI.MEMORY_LIMIT
    else
        return MOI.OTHER_ERROR
    end
end

"""
    MOISolutionStatus(st::SolutionStatus)

Convert a Tulip `SolutionStatus` into a `MOI.ResultStatusCode`.
"""
function MOISolutionStatus(st::SolutionStatus)::MOI.ResultStatusCode
    if st == Sln_Unknown
        return MOI.UNKNOWN_RESULT_STATUS
    elseif st == Sln_Optimal || st == Sln_FeasiblePoint
        return MOI.FEASIBLE_POINT
    elseif st == Sln_InfeasiblePoint
        return MOI.INFEASIBLE_POINT
    elseif st == Sln_InfeasibilityCertificate
        return MOI.INFEASIBILITY_CERTIFICATE
    else
        return MOI.OTHER_RESULT_STATUS
    end
end

"""
    _bounds(s)

"""
_bounds(s::MOI.EqualTo{T}) where{T} = s.value, s.value
_bounds(s::MOI.LessThan{T}) where{T}  = T(-Inf), s.upper
_bounds(s::MOI.GreaterThan{T}) where{T}  = s.lower, T(Inf)
_bounds(s::MOI.Interval{T}) where{T}  = s.lower, s.upper

const SCALAR_SETS{T} = Union{
    MOI.LessThan{T},
    MOI.GreaterThan{T},
    MOI.EqualTo{T},
    MOI.Interval{T}
} where{T}

@enum(ObjType, _SINGLE_VARIABLE, _SCALAR_AFFINE)


# ==============================================================================
# ==============================================================================
#
#               S U P P O R T E D    M O I    F E A T U R E S
#
# ==============================================================================
# ==============================================================================

"""
    Optimizer{T}

Wrapper for MOI.
"""
mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    inner::Model{T}

    objective_sense::Union{Nothing,MOI.OptimizationSense}
    _obj_type::Union{Nothing,ObjType}

    # Map MOI Variable/Constraint indices to internal indices
    var_counter::Int  # Should never be reset
    con_counter::Int  # Should never be reset
    var_indices_moi::Vector{MOI.VariableIndex}
    var_indices::Dict{MOI.VariableIndex, Int}
    con_indices_moi::Vector{MOI.ConstraintIndex}
    con_indices::Dict{MOI.ConstraintIndex{MOI.ScalarAffineFunction{T}, <:SCALAR_SETS{T}}, Int}

    # Variable and constraint names
    name2var::Dict{String, Set{MOI.VariableIndex}}
    name2con::Dict{String, Set{MOI.ConstraintIndex}}

    # Keep track of bound constraints
    var2bndtype::Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}

    # Tulip.Model does not record solution time...
    solve_time::Float64

    function Optimizer{T}(;kwargs...) where{T}
        m = new{T}(
            Model{T}(),
            nothing,  # objective_sense
            nothing,  # _obj_type
            # Variable and constraint counters
            0, 0,
            # Index mapping
            MOI.VariableIndex[], Dict{MOI.VariableIndex, Int}(),
            MOI.ConstraintIndex[], Dict{MOI.ConstraintIndex{MOI.ScalarAffineFunction, <:SCALAR_SETS{T}}, Int}(),
            # Name -> index mapping
            Dict{String, Set{MOI.VariableIndex}}(),
            Dict{String, Set{MOI.ConstraintIndex}}(),
            Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}(),
            0.0
        )

        for (k, v) in kwargs
            set_parameter(m.inner, string(k), v)
        end

        return m
    end
end

Optimizer(;kwargs...) = Optimizer{Float64}(;kwargs...)

function MOI.empty!(m::Optimizer)
    # Inner model
    empty!(m.inner)
    m.objective_sense = nothing
    m._obj_type = nothing
    # Reset index mappings
    m.var_indices_moi = MOI.VariableIndex[]
    m.con_indices_moi = MOI.ConstraintIndex[]
    m.var_indices = Dict{MOI.VariableIndex, Int}()
    m.con_indices = Dict{MOI.ConstraintIndex, Int}()

    # Reset name mappings
    m.name2var = Dict{String, Set{MOI.VariableIndex}}()
    m.name2con = Dict{String, Set{MOI.ConstraintIndex}}()

    # Reset bound tracking
    m.var2bndtype  = Dict{MOI.VariableIndex, Set{MOI.ConstraintIndex}}()

    m.solve_time = 0.0
    return nothing
end

function MOI.is_empty(m::Optimizer)
    m.objective_sense === nothing || return false
    m._obj_type === nothing || return false
    m.inner.pbdata.nvar == 0 || return false
    m.inner.pbdata.ncon == 0 || return false

    length(m.var_indices) == 0 || return false
    length(m.var_indices_moi) == 0 || return false
    length(m.con_indices) == 0 || return false
    length(m.con_indices_moi) == 0 || return false

    length(m.name2var) == 0 || return false
    length(m.name2con) == 0 || return false

    length(m.var2bndtype) == 0 || return false

    return true
end

function MOI.optimize!(m::Optimizer)
    t_solve = @elapsed optimize!(m.inner)
    m.solve_time = t_solve
    return nothing
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(dest, src)
end


# ==============================================================================
#           I. Optimizer attributes
# ==============================================================================
# ==============================================================================
#           II. Model attributes
# ==============================================================================
include("./attributes.jl")

# ==============================================================================
#           III. Variables
# ==============================================================================
include("./variables.jl")

# ==============================================================================
#           IV. Constraints
# ==============================================================================
include("./constraints.jl")

# ==============================================================================
#           V. Objective
# ==============================================================================
include("./objective.jl")
