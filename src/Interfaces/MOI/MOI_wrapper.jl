import MathOptInterface
const MOI = MathOptInterface

# ==============================================================================
#           HELPER FUNCTIONS
# ==============================================================================

"""
    MOITerminationStatus(st::TerminationStatus)

Convert a Tulip `TerminationStatus` into a `MOI.TerminationStatusCode`.
"""
function MOITerminationStatus(st::TerminationStatus)::MOI.TerminationStatusCode
    if st == Trm_Optimal
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
_bounds(s::MOI.EqualTo{Tv}) where{Tv} = s.value, s.value
_bounds(s::MOI.LessThan{Tv}) where{Tv}  = Tv(-Inf), s.upper
_bounds(s::MOI.GreaterThan{Tv}) where{Tv}  = s.lower, Tv(Inf)
_bounds(s::MOI.Interval{Tv}) where{Tv}  = s.lower, s.upper

const SCALAR_SETS{Tv} = Union{
    MOI.LessThan{Tv},
    MOI.GreaterThan{Tv},
    MOI.EqualTo{Tv},
    MOI.Interval{Tv}
} where{Tv}



# ==============================================================================
# ==============================================================================
#
#               S U P P O R T E D    M O I    F E A T U R E S
#
# ==============================================================================
# ==============================================================================

"""
    Optimizer{Tv}

Wrapper for MOI.
"""
mutable struct Optimizer{Tv} <: MOI.AbstractOptimizer
    inner::Model{Tv}

    is_feas::Bool  # Model is feasibility problem if true

    # Map MOI Variable/Constraint indices to internal indices
    var_counter::Int  # Should never be reset
    con_counter::Int  # Should never be reset
    var_indices_moi::Vector{MOI.VariableIndex}
    var_indices::Dict{MOI.VariableIndex, Int}
    con_indices_moi::Vector{MOI.ConstraintIndex}
    con_indices::Dict{MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, <:SCALAR_SETS{Tv}}, Int}

    # Variable and constraint names
    name2var::Dict{String, MOI.VariableIndex}
    name2con::Dict{String, MOI.ConstraintIndex}
    # MOIIndex -> name mapping for SingleVariable constraints
    # Will be dropped with MOI 0.10
    #   => (https://github.com/JuliaOpt/MathOptInterface.jl/issues/832)
    bnd2name::Dict{MOI.ConstraintIndex, String}

    # Keep track of bound constraints
    var2bndtype::Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}

    function Optimizer{Tv}(;kwargs...) where{Tv}
        m = new{Tv}(
            Model{Tv}(), false,
            # Variable and constraint counters
            0, 0, 
            # Index mapping
            MOI.VariableIndex[], Dict{MOI.VariableIndex, Int}(),
            MOI.ConstraintIndex[], Dict{MOI.ConstraintIndex{MOI.ScalarAffineFunction, <:SCALAR_SETS{Tv}}, Int}(),
            # Name -> index mapping
            Dict{String, MOI.VariableIndex}(), Dict{String, MOI.ConstraintIndex}(),
            Dict{MOI.ConstraintIndex, String}(),  # Variable bounds tracking
            Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}()
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
    # Reset index mappings
    m.var_indices_moi = MOI.VariableIndex[]
    m.con_indices_moi = MOI.ConstraintIndex[]
    m.var_indices = Dict{MOI.VariableIndex, Int}()
    m.con_indices = Dict{MOI.ConstraintIndex, Int}()

    # Reset name mappings
    m.name2var = Dict{String, MOI.VariableIndex}()
    m.name2con = Dict{String, MOI.ConstraintIndex}()

    # Reset bound tracking
    m.bnd2name = Dict{MOI.ConstraintIndex, String}()
    m.var2bndtype  = Dict{MOI.VariableIndex, Set{MOI.ConstraintIndex}}()
end

function MOI.is_empty(m::Optimizer)
    m.inner.pbdata.nvar == 0 || return false
    m.inner.pbdata.ncon == 0 || return false

    length(m.var_indices) == 0 || return false
    length(m.var_indices_moi) == 0 || return false
    length(m.con_indices) == 0 || return false
    length(m.con_indices_moi) == 0 || return false

    length(m.name2var) == 0 || return false
    length(m.name2con) == 0 || return false

    length(m.bnd2name) == 0 || return false
    length(m.var2bndtype) == 0 || return false

    return true
end

MOI.optimize!(m::Optimizer) = optimize!(m.inner)

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
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