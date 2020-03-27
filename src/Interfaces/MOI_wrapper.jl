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
_bounds(s::MOI.EqualTo{Tv}) where{Tv<:Real} = s.value, s.value
_bounds(s::MOI.LessThan{Tv}) where{Tv<:Real}  = Tv(-Inf), s.upper
_bounds(s::MOI.GreaterThan{Tv}) where{Tv<:Real}  = s.lower, Tv(Inf)
_bounds(s::MOI.Interval{Tv}) where{Tv<:Real}  = s.lower, s.upper

const SCALAR_SETS{Tv} = Union{
    MOI.LessThan{Tv},
    MOI.GreaterThan{Tv},
    MOI.EqualTo{Tv},
    MOI.Interval{Tv}
} where{Tv<:Real}



# ==============================================================================
# ==============================================================================
#
#               S U P P O R T E D    M O I    F E A T U R E S
#
# ==============================================================================
# ==============================================================================

"""
    Optimizer{Tv<:Real}

Wrapper for MOI.
"""
mutable struct Optimizer{Tv<:Real} <: MOI.AbstractOptimizer
    inner::Model{Tv}

    is_feas::Bool  # Model is feasibility problem if true

    # Keep track of names of variable bound constraints
    # Should go away in future MOI release
    bnd2name::Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}
    name2bnd_LT::Dict{String, MOI.VariableIndex}
    name2bnd_GT::Dict{String, MOI.VariableIndex}
    name2bnd_ET::Dict{String, MOI.VariableIndex}
    name2bnd_IT::Dict{String, MOI.VariableIndex}

    # Keep track of bound constraints
    var2bndtype::Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}

    function Optimizer{Tv}() where{Tv<:Real}
        return new{Tv}(
            Model{Tv}(), false,
            Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{MOI.VariableIndex, Set{Type{<:MOI.AbstractScalarSet}}}()
        )
    end
end

# TODO: set parameters via kw arguments
Optimizer() = Optimizer{Float64}()

function MOI.empty!(m::Optimizer{Tv}) where{Tv<:Real}
    empty!(m.inner)
    m.bnd2name = Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}()
    m.name2bnd_LT = Dict{String, MOI.VariableIndex}()
    m.name2bnd_GT = Dict{String, MOI.VariableIndex}()
    m.name2bnd_ET = Dict{String, MOI.VariableIndex}()
    m.name2bnd_IT = Dict{String, MOI.VariableIndex}()
    m.var2bndtype  = Dict{MOI.VariableIndex, Set{MOI.ConstraintIndex}}()
end

MOI.is_empty(m::Optimizer) = is_empty(m.inner)

MOI.optimize!(m::Optimizer) = optimize!(m.inner)

MOI.Utilities.supports_default_copy_to(::Optimizer, ::Bool) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
end


# ==============================================================================
#           I. Optimizer attributes
# ==============================================================================

    # =============================================
    #   I.1 Supported attributes
    # =============================================
const SUPPORTED_OPTIMIZER_ATTR = Union{
    MOI.NumberOfThreads,
    MOI.SolverName,
    MOI.Silent,
    MOI.TimeLimitSec,
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_OPTIMIZER_ATTR} = true

    # =============================================
    #   I.2 Get/set optimizer attributes
    # =============================================

# =============================================
#   NumberOfThreads
# =============================================
MOI.get(m::Optimizer, ::MOI.NumberOfThreads) = m.inner.env.threads

function MOI.set(m::Optimizer, ::MOI.NumberOfThreads, n::Int)
    m.inner.env.threads = n
    return nothing
end

# =============================================
#   SolverName
# =============================================
MOI.get(::Optimizer, ::MOI.SolverName) = "Tulip"

# =============================================
#   Silent
# =============================================
MOI.get(m::Optimizer, ::MOI.Silent) = iszero(m.inner.env.verbose)

function MOI.set(m::Optimizer, ::MOI.Silent, flag::Bool)
    m.inner.env.verbose = 1 - flag
    return nothing
end

# =============================================
#   TimeLimitSec
# =============================================
MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m.inner.env.time_limit

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, t)
    m.inner.env.time_limit = t
    return nothing
end

# ==============================================================================
#           II. Model attributes
# ==============================================================================

    # =============================================
    #   II.1 Supported attributes
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
    # MOI.ObjectiveBound,  # This is generally used for MIP models
    MOI.RelativeGap,
    # MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    # MOI.RawStatusString,  # TODO
    MOI.ResultCount,
    MOI.TerminationStatus,
    MOI.PrimalStatus,
    MOI.DualStatus
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_MODEL_ATTR} = true

# =============================================
#   Name
# =============================================
MOI.get(m::Optimizer, ::MOI.Name) = m.inner.name

# =============================================
#   ObjectiveSense
# =============================================
function MOI.get(m::Optimizer, ::MOI.ObjectiveSense)
    m.is_feas && return MOI.FEASIBILITY_SENSE

    s = m.inner.pbdata_raw.obj_sense
    if s == TLP_MIN
        return MOI.MIN_SENSE
    else
        return MOI.MAX_SENSE
    end
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense)
    m.is_feas = (s == MOI.FEASIBILITY_SENSE)

    if s == MOI.MIN_SENSE || s == MOI.FEASIBILITY_SENSE
        m.inner.pbdata_raw.obj_sense = TLP_MIN
    elseif s == MOI.MAX_SENSE
        m.inner.pbdata_raw.obj_sense = TLP_MAX
    else
        error("Objetive sense not supported: $s")
    end
    
    return nothing
end

# =============================================
#   NumberOfVariable
# =============================================
MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = get_num_var(m.inner)

# =============================================
#   ListOfVariableIndices
# =============================================
function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return [
        MOI.VariableIndex(vidx.uuid)
        for vidx in keys(m.inner.pbdata_raw.vars)
    ]
end

# =============================================
#   ListOfConstraintIndices
# =============================================
function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ListOfConstraintIndices{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[]

    var_indices = MOI.get(m, MOI.ListOfVariableIndices())

    for var in var_indices
        if S ∈ m.var2bndtype[var]
            push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(var.value))
        end
    end
    return indices
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ListOfConstraintIndices{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    indices = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}[]

    for (cidx, con) in m.inner.pbdata_raw.constrs
        bt, lb, ub = get_bounds(con)
        if (
            (S == MOI.LessThan{Tv} && (bt == TLP_UP))
            || (S == MOI.GreaterThan{Tv} && (bt == TLP_LO))
            || (S == MOI.EqualTo{Tv} && (bt == TLP_FX))
            || (S == MOI.Interval{Tv} && (bt == TLP_RG))
        )
            push!(indices, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid))
        end
    end
    return indices
end

# =============================================
#   NumberOfConstraints
# =============================================
function MOI.get(
    m::Optimizer{Tv},
    ::MOI.NumberOfConstraints{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    ncon = 0

    for v in MOI.get(m, MOI.ListOfVariableIndices())
        ncon += S ∈ m.var2bndtype[v]
    end

    return ncon
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    ncon = 0

    for (cidx, con) in m.inner.pbdata_raw.constrs
        bt, lb, ub = get_bounds(con)
        if (
            (S == MOI.LessThan{Tv} && (bt == TLP_UP))
            || (S == MOI.GreaterThan{Tv} && (bt == TLP_LO))
            || (S == MOI.EqualTo{Tv} && (bt == TLP_FX))
            || (S == MOI.Interval{Tv} && (bt == TLP_RG))
        )
            ncon += 1
        end
    end
    return ncon
end

# =============================================
#   ObjectiveFunctionType
# =============================================
MOI.get(
    m::Optimizer{Tv}, ::MOI.ObjectiveFunctionType
) where{Tv<:Real} = MOI.ScalarAffineFunction{Tv}

# =============================================
#   ObjectiveValue
# =============================================
function MOI.get(m::Optimizer{Tv}, attr::MOI.ObjectiveValue) where{Tv<:Real}
    MOI.check_result_index_bounds(m, attr)

    # TODO: dispatch a function call on m.inner
    if m.inner.pbdata_raw.obj_sense == TLP_MIN
        return m.inner.solver.primal_bound_scaled
    elseif m.inner.pbdata_raw.obj_sense == TLP_MAX
        return - m.inner.solver.primal_bound_scaled
    end
end

# =============================================
#   DualObjectiveValue
# =============================================
function MOI.get(m::Optimizer{Tv}, attr::MOI.DualObjectiveValue) where{Tv<:Real}
    MOI.check_result_index_bounds(m, attr)
    
    # TODO: dispatch a function call on m.inner
    if m.inner.pbdata_raw.obj_sense == TLP_MIN
        return m.inner.solver.dual_bound_scaled
    elseif m.inner.pbdata_raw.obj_sense == TLP_MAX
        return - m.inner.solver.primal_bound_scaled
    end
end

# =============================================
#   RelativeGap
# =============================================
function MOI.get(m::Optimizer{Tv}, ::MOI.RelativeGap) where{Tv<:Real}
    # TODO: dispatch a function call on m.inner
    zp = m.inner.solver.primal_bound_scaled
    zd = m.inner.solver.dual_bound_scaled
    return (abs(zp - zd) / (Tv(1e-10) + abs(zd)))
end

# =============================================
#   SolveTime
# =============================================
# TODO: MOI.get(m::Optimizer, ::MOI.SolveTime)

# =============================================
#   SimplexIterations
# =============================================
MOI.get(::Optimizer, ::MOI.SimplexIterations) = 0

# =============================================
#   BarrierIterations
# =============================================
# TODO: dispatch function call on m.inner
MOI.get(m::Optimizer, ::MOI.BarrierIterations) = m.inner.solver.niter

# =============================================
#   RawSolver
# =============================================
MOI.get(m::Optimizer, ::MOI.RawSolver) = m.inner

# =============================================
#   ResultCount
# =============================================
function MOI.get(m::Optimizer, ::MOI.ResultCount)
    st = MOI.get(m, MOI.TerminationStatus())

    if (
        st == MOI.OPTIMIZE_NOT_CALLED
        || st == MOI.OTHER_ERROR
        || st == MOI.MEMORY_LIMIT
    )
        return 0
    end
    return 1
end

# =============================================
#   TerminationStatus
# =============================================
function MOI.get(m::Optimizer, ::MOI.TerminationStatus)
    # TODO: dispatch function call on m.inner
    isa(m.inner.solver, Nothing) && return MOI.OPTIMIZE_NOT_CALLED

    st = m.inner.solver.solver_status
    
    return MOITerminationStatus(st)
end

# =============================================
#   PrimalStatus
# =============================================
function MOI.get(m::Optimizer, attr::MOI.PrimalStatus)
    attr.N == 1 || return MOI.NO_SOLUTION

    return MOISolutionStatus(m.inner.solver.primal_status)
end

# =============================================
#   DualStatus
# =============================================
function MOI.get(m::Optimizer, attr::MOI.DualStatus)
    attr.N == 1 || return MOI.NO_SOLUTION
    
    return MOISolutionStatus(m.inner.solver.dual_status)
end


# ==============================================================================
#           III. Variables
# ==============================================================================

    # =============================================
    #   III.1 Supported variable attributes
    # =============================================
"""
    SUPPORTED_VARIABLE_ATTR

List of supported MOI `VariableAttribute`.
* `MOI.VariableName`
* `MOI.VariablePrimal`
"""
const SUPPORTED_VARIABLE_ATTR = Union{
    MOI.VariableName,
    # MOI.VariablePrimalStart,
    MOI.VariablePrimal
}

MOI.supports(::Optimizer, ::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

    # =============================================
    #   III.2 Add variables
    # =============================================
function MOI.add_variable(m::Optimizer{Tv}) where{Tv<:Real}
    # By default, we create a variable as follows:
    # * empty name
    # * zero objective coeff
    # * lower and upper bounds are -/+ infinity (i.e., free variable)
    # * does not appear in any constraint
    # TODO: dispatch a function call to m.inner
    vidx = Tulip.add_variable!(m.inner, "", zero(Tv), Tv(-Inf), Tv(Inf), ConstrId[], Tv[])
    m.var2bndtype[MOI.VariableIndex(vidx.uuid)] = Set{MOI.ConstraintIndex}()
    return MOI.VariableIndex(vidx.uuid)
end

# TODO: dispatch to inner model
function MOI.add_variables(m::Optimizer{Tv}, N::Int) where{Tv<:Real}
    N >= 0 || error("Cannot add negative number of variables")

    vars = MOI.VariableIndex[]
    for i in 1:N
        x = MOI.add_variable(m)
        push!(vars, x)
    end

    return vars
end

function MOI.is_valid(m::Optimizer, v::MOI.VariableIndex)
    return haskey(m.inner.pbdata_raw.vars, VarId(v.value))
end

    # =============================================
    #   III.3 Delete variables
    # =============================================
function MOI.delete(m::Optimizer{Tv}, v::MOI.VariableIndex) where{Tv<:Real}
    MOI.throw_if_not_valid(m, v)

    # Remove variable from model
    delete_variable!(m.inner.pbdata_raw, VarId(v.value))
    
    # Remove bound tracking
    delete!(m.var2bndtype, v)
end

    # =============================================
    #   III.4 Get/set variable attributes
    # =============================================

# =============================================
#   VariableIndex
# =============================================
function MOI.get(
    m::Optimizer{Tv},
    ::Type{MOI.VariableIndex},
    name::String
) where{Tv<:Real}
    if haskey(m.inner.pbdata_raw.name2var, name)
        vidx = m.inner.pbdata_raw.name2var[name]
        return MOI.VariableIndex(vidx.uuid)
    else    
        # TODO: returning nothing is not type-stable.
        # Raise the issue in MOI
        return nothing
    end
end

# =============================================
#   VariableName
# =============================================
function MOI.get(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex)
    MOI.is_valid(m, v) || throw(MOI.InvalidIndex(v))
    return get_var_name(m.inner, VarId(v.value))
end

function MOI.set(m::Optimizer, ::MOI.VariableName, v::MOI.VariableIndex, name::String)
    # Check that variable does exist
    MOI.is_valid(m, v) || throw(MOI.InvalidIndex(v))

    # Check that name is unique
    !haskey(m.inner.pbdata_raw.name2var, name) || error("Variable $name already exists.")

    # TODO: add a function at the ProblemData level
    vidx = VarId(v.value)
    var = m.inner.pbdata_raw.vars[vidx]
    set_name!(var, name)
    # Update name in dict
    m.inner.pbdata_raw.name2var[name] = vidx
    return nothing
end

# =============================================
#   VariablePrimal
# =============================================
function MOI.get(
    m::Optimizer{Tv}, attr::MOI.VariablePrimal,
    x::MOI.VariableIndex
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, x)

    # Check result index
    MOI.check_result_index_bounds(m, attr)

    return get_value(m.inner, VarId(x.value))
end


# ==============================================================================
#           IV. Constraints
# ==============================================================================

    # =============================================
    #   IV.1 Supported attributes and constraints
    # =============================================

"""
    SUPPORTED_CONSTR_ATTR

List of supported MOI `ConstraintAttribute`.
* `MOI.ConstraintName`
"""
const SUPPORTED_CONSTR_ATTR = Union{
    MOI.ConstraintName,
    MOI.ConstraintPrimal,
    MOI.ConstraintDual,
    # MOI.ConstraintPrimalStart,
    # MOI.ConstraintDualStart, # once dual warm-start is supported 
    # MOI.ConstraintBasisStatus,  # once cross-over is supported
    MOI.ConstraintFunction,
    MOI.ConstraintSet
}

MOI.supports(::Optimizer, ::A, ::Type{<:MOI.ConstraintIndex}) where{A<:SUPPORTED_CONSTR_ATTR} = true

# Variable bounds
# Hopefully, MOI conventions will change and we can get rid of all the extra code
# for SingleVariable-In-LessThan/GreaterThan/EqualTo/Interval constraints
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{S}
) where {Tv<:Real, S<:SCALAR_SETS{Tv}}
    return true
end

# General linear constraints
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Tv}}, ::Type{S}
) where {Tv<:Real, S<:SCALAR_SETS{Tv}}
    return true
end

    # =============================================
    #   IV.2 Add constraints
    # =============================================

# TODO: make it clear that only finite bounds can be given in input.
# To relax variable bounds, one needs to delete the associated bound constraint.
function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.LessThan{Tv}
) where{Tv<:Real}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    
    # identify which variable
    var = m.inner.pbdata_raw.vars[VarId(v.value)]

    # Check if upper bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.LessThan{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.LessThan{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.LessThan{Tv}}(v.value))
    end

    # Get current bounds
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}(v.value)
    bt, lb, ub = get_bounds(var)

    # Update upper-bound
    set_bounds!(var, lb, s.upper)
    push!(m.var2bndtype[v], MOI.LessThan{Tv})
    
    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.GreaterThan{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    
    # identify which variable
    var = m.inner.pbdata_raw.vars[VarId(v.value)]

    # Check if lower bound already exists
    if MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.GreaterThan{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Tv}, MOI.GreaterThan{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Tv}, MOI.GreaterThan{Tv}}(v.value))
    end

    # Get current bounds
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}(v.value)
    bt, lb, ub = get_bounds(var)

    # Update upper-bound
    set_bounds!(var, s.lower, ub)
    push!(m.var2bndtype[v], MOI.GreaterThan{Tv})

    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.EqualTo{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    
    # identify which variable
    var = m.inner.pbdata_raw.vars[VarId(v.value)]

    # Check if a bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.EqualTo{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.EqualTo{Tv}}(v.value))
    end

    # Update variable bounds
    set_bounds!(var, s.value, s.value)
    push!(m.var2bndtype[v], MOI.EqualTo{Tv})

    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}(v.value)

    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.Interval{T}
) where{T<:Real, Tv<:Real}
    # Check that variable exists
    v = f.variable
    MOI.throw_if_not_valid(m, v)
    
    # identify which variable
    var = m.inner.pbdata_raw.vars[VarId(v.value)]

    # Check if a bound already exists
    if MOI.LessThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.GreaterThan{Tv} ∈ m.var2bndtype[v]
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.EqualTo{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.EqualTo{Tv}, MOI.Interval{Tv}}(v.value))
    elseif MOI.Interval{Tv} ∈ m.var2bndtype[v]
        throw(MOI.UpperBoundAlreadySet{MOI.Interval{Tv}, MOI.Interval{Tv}}(v.value))
    end

    # Update variable bounds
    set_bounds!(var, s.lower, s.upper)
    push!(m.var2bndtype[v], MOI.Interval{Tv})

    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}(v.value)
    return cidx
end

# General linear constraints
function MOI.add_constraint(
    m::Optimizer,
    f::MOI.ScalarAffineFunction{Tv},
    s::S
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    # Check that constant term is zero
    if !iszero(f.constant)
        throw(MOI.ScalarFunctionConstantNotZero{Tv, typeof(f), S}(f.constant))
    end

    # Extract bound info
    # TODO: check that bounds are finite
    lb, ub = _bounds(s)

    # Convert to canonical form
    f_canon = MOI.Utilities.canonical(f)
    # Create new constraint
    colids = [
        VarId(t.variable_index.value)
        for t in f_canon.terms
    ]
    colvals = [
        t.coefficient
        for t in f_canon.terms
    ]

    cidx = add_constraint!(m.inner, "", lb, ub, colids, colvals)

    # Return MOI index
    return MOI.ConstraintIndex{typeof(f), S}(cidx.uuid)
end

    # =============================================
    #   IV.3 Index checking
    # =============================================
function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S <:SCALAR_SETS{Tv}}
    
    v = MOI.VariableIndex(c.value)
    MOI.is_valid(m, v) || return false
    
    return S ∈ m.var2bndtype[v]
end

function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    # Get constraint index
    cidx = ConstrId(c.value)

    # Check if constraint exists in model
    if !haskey(m.inner.pbdata_raw.constrs, cidx)
        # Constraint does not exist in model
        return false
    end

    # Check that bound type matches that of `c`
    bt, lb, ub = get_constr_bounds(m.inner, cidx)
    return (
           (S == MOI.LessThan{Tv}       && bt == TLP_UP)
        || (S == MOI.GreaterThan{Tv}    && bt == TLP_LO)
        || (S == MOI.EqualTo{Tv}        && bt == TLP_FX)
        || (S == MOI.Interval{Tv}       && bt == TLP_RG)
    )
end

    # =============================================
    #   IV.4 Delete constraints
    # =============================================  
function MOI.delete(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    v = MOI.VariableIndex(c.value)

    # Check if constraint had a name. If so, it will need to be deleted.
    name = MOI.get(m, MOI.ConstraintName(), c)

    # Get bounds
    var = m.inner.pbdata_raw.vars[VarId(c.value)]
    bt, lb, ub = get_bounds(var)

    # Update bounds
    if S == MOI.LessThan{Tv}
        # Remove upper-bound
        set_bounds!(var, lb, Tv(Inf))
    elseif S == MOI.GreaterThan{Tv}
        # Remove lower bound
        set_bounds!(var, Tv(-Inf), ub)
    else
        # Set variable to free
        set_bounds!(var, Tv(-Inf), Tv(Inf))
    end

    # Delete bound constraint name
    delete!(m.bnd2name, c)
    if name != ""
        if S == MOI.LessThan{Tv}
            delete!(m.name2bnd_LT, v)

        elseif S == MOI.GreaterThan{Tv}
            delete!(m.name2bnd_GT, v)

        elseif S == MOI.EqualTo{Tv}
            delete!(m.name2bnd_ET, v)

        elseif S == MOI.Interval{Tv}
            delete!(m.name2bnd_IT, v)
        end

    end

    # Delete tracking of bounds
    delete!(m.var2bndtype[v], S)

    return nothing
end

function MOI.delete(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    cidx = ConstrId(c.value)

    delete_constraint!(m.inner.pbdata_raw, cidx)
end

    # =============================================
    #   IV.5 Modify constraints
    # =============================================
function MOI.modify(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    chg::MOI.ScalarCoefficientChange{Tv}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))
    MOI.is_valid(m, chg.variable) || throw(MOI.InvalidIndex(chg.variable))

    set_coeff!(
        m.inner.pbdata_raw,
        VarId(chg.variable.value),
        ConstrId(c.value),
        chg.new_coefficient
    )
    return nothing
end

    # =============================================
    #   IV.4 Get/set constraint attributes
    # ============================================= 

# =============================================
#   ConstraintName
# =============================================
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    
    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Get constraint name for dict directly
    if !haskey(m.bnd2name, c)
        return ""
    else
        return m.bnd2name[c]
    end
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    
    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Get constraint index
    cidx = ConstrId(c.value)
    return get_constr_name(m.inner, cidx)
end

# Set name from constraint index
function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    name::String
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}

    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Set constraint name in dict directly
    m.bnd2name[c] = name
    if S == MOI.LessThan{Tv}
        m.name2bnd_LT[name] = MOI.VariableIndex(c.value)

    elseif S == MOI.GreaterThan{Tv}
        m.name2bnd_GT[name] = MOI.VariableIndex(c.value)

    elseif S == MOI.EqualTo{Tv}
        m.name2bnd_ET[name] = MOI.VariableIndex(c.value)

    elseif S == MOI.Interval{Tv}
        m.name2bnd_IT[name] = MOI.VariableIndex(c.value)
    end
    return nothing
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    name::String
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}

    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Get constraint index
    cidx = ConstrId(c.value)
    con = m.inner.pbdata_raw.constrs[cidx]
    con.dat.name = name

    # Update name in dict
    m.inner.pbdata_raw.name2con[name] = cidx
    return nothing
end

# =============================================
#   ConstraintIndex
# =============================================
# Get constraint index from name
function MOI.get(
    m::Optimizer{Tv},
    ::Type{MOI.ConstraintIndex},
    name::String
) where{Tv<:Real}
    for S in [MOI.LessThan{Tv}, MOI.GreaterThan{Tv}, MOI.EqualTo{Tv}, MOI.Interval{Tv}]
        # Check in linear constraints.
        cidx = MOI.get(m, MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}, name)
        isa(cidx, Nothing) || return cidx

        # Check in bound constraints
        cidx = MOI.get(m, MOI.ConstraintIndex{MOI.SingleVariable, S}, name)
        isa(cidx, Nothing) || return cidx
    end

    # Could not find name anywhere
    return nothing
end

function MOI.get(
    m::Optimizer{Tv},
    ::Type{MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}},
    name::String
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    if haskey(m.inner.pbdata_raw.name2con, name)
        # Extract MOI index and return
        cidx = m.inner.pbdata_raw.name2con[name]

        # check bound
        bt, lb, ub = get_constr_bounds(m.inner, cidx)

        if bt == TLP_UP && S == MOI.LessThan{Tv}
            return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid)
        elseif bt == TLP_LO && S == MOI.GreaterThan{Tv}
            return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid)
        elseif bt == TLP_FX && S == MOI.EqualTo{Tv}
            return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid)
        elseif bt == TLP_RG && S == MOI.Interval{Tv}
            return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid)
        end
    else
        return nothing
    end
end

function MOI.get(
    m::Optimizer{Tv},
    ::Type{MOI.ConstraintIndex{MOI.SingleVariable, S}},
    name::String
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    # check if name exists
    if S == MOI.LessThan{Tv}
        haskey(m.name2bnd_LT, name) || return nothing
        vid = m.name2bnd_LT[name]

    elseif S == MOI.GreaterThan{Tv}
        haskey(m.name2bnd_GT, name) || return nothing
        vid = m.name2bnd_GT[name]

    elseif S == MOI.EqualTo{Tv}
        haskey(m.name2bnd_ET, name) || return nothing
        vid = m.name2bnd_ET[name]

    elseif S == MOI.Interval{Tv}
        haskey(m.name2bnd_IT, name) || return nothing
        vid = m.name2bnd_IT[name]

    else
        return nothing
    end

    return MOI.ConstraintIndex{MOI.SingleVariable, S}(vid.value)

end
    
# =============================================
#   ConstraintPrimal
# ============================================= 
function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    # Check result index
    MOI.check_result_index_bounds(m, attr)

    return MOI.get(m, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)
    MOI.check_result_index_bounds(m, attr)

    return MOI.Utilities.get_fallback(m, MOI.ConstraintPrimal(), c)
end

# =============================================
#   ConstraintDual
# =============================================
function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check
    MOI.check_result_index_bounds(m, attr)

    # Get variable index
    vidx = VarId(c.value)
    v = MOI.VariableIndex(c.value)

    # Get variable bounds
    bt, lb, ub = get_var_bounds(m.inner, vidx)
    # Get index in standard form
    vind = m.inner.pbdata_std.var2idx[vidx]

    if (MOI.LessThan{Tv} ∈ m.var2bndtype[v]) && (MOI.GreaterThan{Tv} ∉ m.var2bndtype[v])
        # x_i <= u was transformed into -x_i >= -u
        # Therefore, we need to return -s_i

        # Get corresponding dual
        s = m.inner.solver.pt.s[vind] / m.inner.solver.pt.t
        return -s
    elseif (MOI.LessThan{Tv} ∈ m.var2bndtype[v]) && (MOI.GreaterThan{Tv} ∈ m.var2bndtype[v])
        # l <= x <= u
        # We need to return the index that corresponds to the upper-bound constraint

        # Get index of upper-bound constraint
        # TODO: have a more convenient way of doing this
        for (i, j) in enumerate(m.inner.pbdata_std.uind)
            if j == vind
                # bingo
                s = m.inner.solver.pt.z[i] / m.inner.solver.pt.t
                return -s
            end
        end
    end
    error("Could not find dual for upper-bound on variable $vind.")
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check
    MOI.check_result_index_bounds(m, attr)

    # Get variable index
    vidx = VarId(c.value)
    
    # Get index in standard form
    vind = m.inner.pbdata_std.var2idx[vidx]

    # Get corresponding dual
    s = m.inner.solver.pt.s[vind] / m.inner.solver.pt.t
    return s
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:Union{MOI.EqualTo{Tv}, MOI.Interval{Tv}}}
    MOI.throw_if_not_valid(m, c)  # Sanity check
    MOI.check_result_index_bounds(m, attr)

    # Constraint is l <= x <= u
    # The corresponding dual is s_l - s_u

    # Get variable index
    vidx = VarId(c.value)

    # Get index in standard form
    vind = m.inner.pbdata_std.var2idx[vidx]

    # Get s_l
    s_l = m.inner.solver.pt.s[vind] / m.inner.solver.pt.t

    # Get s_u
    s_u = zero(Tv)
    for (i, j) in enumerate(m.inner.pbdata_std.uind)
        if j == vind
            # bingo
            s_u = m.inner.solver.pt.z[i] / m.inner.solver.pt.t
            break
        end
    end

    return s_l - s_u
end

function MOI.get(
    m::Optimizer{Tv}, attr::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check
    MOI.check_result_index_bounds(m, attr)

    return get_dual_value(m.inner, ConstrId(c.value))
end

# =============================================
#   ConstraintFunction
# ============================================= 
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    return MOI.SingleVariable(MOI.VariableIndex(c.value))
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    ::MOI.SingleVariable
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    return throw(MOI.SettingSingleVariableFunctionNotAllowed())
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get constraint index
    cidx = ConstrId(c.value)
    con = m.inner.pbdata_raw.constrs[cidx]

    # Extract constraint coefficients
    terms = MOI.ScalarAffineTerm{Tv}[]
    for vidx in m.inner.pbdata_raw.con2var[cidx]
        # get corresponding coeff
        t = MOI.ScalarAffineTerm(
            m.inner.pbdata_raw.coeffs[vidx, cidx],
            MOI.VariableIndex(vidx.uuid)
        )
        push!(terms, t)
    end

    return MOI.ScalarAffineFunction(terms, zero(Tv))
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintFunction,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    f::MOI.ScalarAffineFunction{Tv}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check
    iszero(f.constant) || throw(
        MOI.ScalarFunctionConstantNotZero{Tv, typeof(f), S}(f.constant)
    )

    # Get constraint index
    cidx = ConstrId(c.value)
    con = m.inner.pbdata_raw.constrs[cidx]

    # Consolidate duplicate terms and remove zeros
    f_canon = MOI.Utilities.canonical(f)

    # Set old row to zero
    f_old = MOI.get(m, MOI.ConstraintFunction(), c)
    for term in f_old.terms
        vidx = VarId(term.variable_index.value)
        set_coeff!(m.inner.pbdata_raw, vidx, cidx, zero(Tv))
    end

    # Set new coeffs
    for term in f_canon.terms
        vidx = VarId(term.variable_index.value)
        val = term.coefficient
        set_coeff!(m.inner.pbdata_raw, vidx, cidx, val)
    end

    # Done

    return nothing
end

# =============================================
#   ConstraintSet
# ============================================= 
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable bounds
    vidx = VarId(c.value)
    bt, lb, ub = get_var_bounds(m.inner, vidx)

    return MOI.LessThan(ub)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable bounds
    vidx = VarId(c.value)
    bt, lb, ub = get_var_bounds(m.inner, vidx)

    return MOI.GreaterThan(lb)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable bounds
    vidx = VarId(c.value)
    bt, lb, ub = get_var_bounds(m.inner, vidx)

    return MOI.EqualTo(lb)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable bounds
    vidx = VarId(c.value)
    bt, lb, ub = get_var_bounds(m.inner, vidx)

    return MOI.Interval(lb, ub)
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get constraint
    cidx = ConstrId(c.value)
    bt, lb, ub = get_constr_bounds(m.inner, cidx)

    if S == MOI.LessThan{Tv}
        return MOI.LessThan(ub)
    elseif S == MOI.GreaterThan{Tv}
        return MOI.GreaterThan(lb)
    elseif S == MOI.EqualTo{Tv}
        return MOI.EqualTo(lb)
    elseif S == MOI.Interval{Tv}
        return MOI.Interval(lb, ub)
    end
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    s::S
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    # Get variable
    vidx = VarId(c.value)
    var = m.inner.pbdata_raw.vars[vidx]

    # Get current bounds
    bt, lb, ub = get_bounds(var)

    # Update upper bound
    # Bound key does not need to be updated
    if S == MOI.LessThan{Tv}
        set_bounds!(var, lb, s.upper)
    elseif S == MOI.GreaterThan{Tv}
        set_bounds!(var, s.lower, ub)
    elseif S == MOI.EqualTo{Tv}
        set_bounds!(var, s.value, s.value)
    elseif S == MOI.Interval{Tv}
        set_bounds!(var, s.lower, s.upper)
    else
        error("Unknown type for ConstraintSet: $S.")
    end

    return nothing
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintSet,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S},
    s::S
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    # Get constraint id
    cid = ConstrId(c.value)
    con = m.inner.pbdata_raw.constrs[cid]

    # Update bounds
    bt, lb, ub = get_bounds(con)

    if S == MOI.LessThan{Tv}
        set_bounds!(con, lb, s.upper)
    elseif S == MOI.GreaterThan{Tv}
        set_bounds!(con, s.lower, ub)
    elseif S == MOI.EqualTo{Tv}
        set_bounds!(con, s.value, s.value)
    elseif S == MOI.Interval{Tv}
        set_bounds!(con, s.lower, s.upper)
    else
        error("Unknown type for ConstraintSet: $S.")
    end

    return nothing
end


# ==============================================================================
#           V. Objective
# ==============================================================================
    
    # =============================================
    #   V.1 Supported objectives
    # =============================================
function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{F}
) where {Tv<:Real, F<:Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Tv}}}
    return true
end

    # =============================================
    #   V.2 Get/set objective function
    # =============================================
function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.SingleVariable}
) where{Tv<:Real}
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}())
    return convert(MOI.SingleVariable, obj)
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}
) where{Tv<:Real}
    
    # Objective coeffs
    terms = MOI.ScalarAffineTerm{Tv}[]
    for (vid, var) in m.inner.pbdata_raw.vars
        c = var.dat.obj
        if !iszero(c)
            push!(terms, MOI.ScalarAffineTerm{Tv}(c, MOI.VariableIndex(vid.uuid)))
        end
    end

    # Constant term
    c0 = m.inner.pbdata_raw.obj_const

    return MOI.ScalarAffineFunction{Tv}(terms, c0)
end

function MOI.get(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Tv}}
) where{Tv<:Real}
    obj = MOI.get(m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}())
    return convert(MOI.ScalarQuadraticFunction{Tv}, obj)
end

function MOI.set(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.SingleVariable},
    f::MOI.SingleVariable
) where{Tv<:Real}

    MOI.set(
        m, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}(),
        convert(MOI.ScalarAffineFunction{Tv}, f)
    )

    return nothing
end

function MOI.set(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    f::MOI.ScalarAffineFunction{T}
) where{T<:Real, Tv<:Real}

    isfinite(f.constant) || error("Objective constant term must be finite.")

    # First, set objective to zero
    for (vidx, var) in m.inner.pbdata_raw.vars
        set_obj_coeff!(var, zero(Tv))
    end

    # Update new coefficients
    # Since there may be dupplicate terms in objective function,
    # we update 
    for t in f.terms
        var = m.inner.pbdata_raw.vars[VarId(t.variable_index.value)]
        c = get_obj_coeff(var)
        set_obj_coeff!(var, c + t.coefficient)
    end

    # Update objective offset
    m.inner.pbdata_raw.obj_const = f.constant

    return nothing
end

    # =============================================
    #   V.3 Modify objective
    # ============================================= 
function MOI.modify(
    m::Optimizer{Tv},
    c::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}},
    chg::MOI.ScalarCoefficientChange{Tv}
) where{Tv<:Real}
    # get index of variable to change
    vidx = VarId(chg.variable.value)

    # Set new coeff
    # TODO: use Tulip API
    var = m.inner.pbdata_raw.vars[vidx]
    set_obj_coeff!(var, chg.new_coefficient)

    return nothing
end

function MOI.modify(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}},
    chg::MOI.ScalarConstantChange{Tv}
) where{Tv<:Real}
    m.inner.pbdata_raw.obj_const = chg.new_constant
    return nothing
end
