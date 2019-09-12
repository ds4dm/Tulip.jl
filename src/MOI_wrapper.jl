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

    # Keep track of names of variable bound constraints
    # Should go away in future MOI release
    bnd2name::Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}
    name2bnd_LT::Dict{String, MOI.VariableIndex}
    name2bnd_GT::Dict{String, MOI.VariableIndex}
    name2bnd_ET::Dict{String, MOI.VariableIndex}
    name2bnd_IT::Dict{String, MOI.VariableIndex}
    function Optimizer{Tv}() where{Tv<:Real}
        return new{Tv}(
            Model{Tv}(),
            Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}(),
            Dict{String, MOI.VariableIndex}()
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
SUPPORTED_OPTIMIZER_ATTR = Union{
    MOI.SolverName,
    MOI.Silent,
    MOI.TimeLimitSec
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_OPTIMIZER_ATTR} = true

    # =============================================
    #   I.2 Get/set optimizer attributes
    # =============================================

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

SUPPORTED_MODEL_ATTR = Union{
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
    MOI.RawStatusString,
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
    s = m.inner.pbdata_raw.obj_sense
    if s == TLP_MIN
        return MOI.MIN_SENSE
    else
        return MOI.MAX_SENSE
    end
end

function MOI.set(m::Optimizer, ::MOI.ObjectiveSense, s::MOI.OptimizationSense)
    
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
        for vidx in keys(m.pbdata_raw.vars)
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

    for (vidx, var) in m.inner.pbdata_raw.vars
        bt, lb, ub = get_bounds(var)
        if (
            (S == MOI.LessThan{Tv} && (bt == TLP_UP || bt == TLP_UL))
            || (S == MOI.GreaterThan{Tv} && (bt == TLP_LO || bt == TLP_UL))
            || (S == MOI.EqualTo{Tv} && (bt == TLP_FX))
            || (S == MOI.Interval{Tv} && (bt == TLP_RG))
        )
            push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(vidx.uuid))
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

    for (vidx, var) in m.inner.pbdata_raw.vars
        bt, lb, ub = get_bounds(var)
        if (
            (S == MOI.LessThan{Tv} && (bt == TLP_UP || bt == TLP_UL))
            || (S == MOI.GreaterThan{Tv} && (bt == TLP_LO || bt == TLP_UL))
            || (S == MOI.EqualTo{Tv} && (bt == TLP_FX))
            || (S == MOI.Interval{Tv} && (bt == TLP_RG))
        )
            ncon += 1
        end
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
function MOI.get(m::Optimizer{Tv}, ::MOI.ObjectiveValue) where{Tv<:Real}
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
function MOI.get(m::Optimizer{Tv}, ::MOI.DualObjectiveValue) where{Tv<:Real}
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
MOI.get(::Optimizer, ::MOI.ResultCount) = 1

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
MOI.get(m::Optimizer, ::MOI.PrimalStatus) = MOISolutionStatus(m.inner.solver.primal_status)

# =============================================
#   DualStatus
# =============================================
MOI.get(m::Optimizer, ::MOI.DualStatus) = MOISolutionStatus(m.inner.solver.dual_status)


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
SUPPORTED_VARIABLE_ATTR = Union{
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

    vidx = VarId(v.value)

    delete_variable!(m.inner.pbdata_raw, vidx)
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
    m::Optimizer{Tv}, ::MOI.VariablePrimal, x::MOI.VariableIndex
) where{Tv<:Real}
    MOI.is_valid(m, x) || throw(MOI.InvalidIndex(x))
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
SUPPORTED_CONSTR_ATTR = Union{
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
    s::MOI.LessThan{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    MOI.is_valid(m, f.variable) || throw(MOI.InvalidIndex(f.variable))

    # identify which variable
    vidx = VarId(f.variable.value)
    var = m.inner.pbdata_raw.vars[vidx]

    # Get current bounds
    bt, lb, ub = get_bounds(var)

    # Check if upper bound already exists
    if bt == TLP_UP || bt == TLP_UL
        S1 = MOI.LessThan{Tv}
        throw(MOI.UpperBoundAlreadySet{S1, typeof(s)}(vidx))
    elseif bt == TLP_FX
        S1 = MOI.EqualTo{Tv}
        throw(MOI.UpperBoundAlreadySet{S1, typeof(s)}(vidx))
    elseif bt == TLP_RG
        S1 = MOI.Interval{Tv}
        throw(MOI.UpperBoundAlreadySet{S1, typeof(s)}(vidx))
    end

    # Upper-bound can be added now
    # Two cases: TLP_UP => TLP_UL; TLP_FR => TLP_UP
    if bt == TLP_LO
        # update bounds
        var.dat.ub = s.upper
        # update bound key
        var.dat.bt = TLP_UL
    elseif bt == TLP_FR
        var.dat.ub = s.upper
        var.dat.bt = TLP_UP
    else
        error("Unsupported bound type: $bt")
    end
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}(vidx.uuid)
    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.GreaterThan{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    MOI.is_valid(m, f.variable) || throw(MOI.InvalidIndex(f.variable))
    
    # identify which variable
    vidx = VarId(f.variable.value)
    var = m.inner.pbdata_raw.vars[vidx]

    # Get current bounds
    bt, lb, ub = get_bounds(var)

    # Check if upper bound already exists
    if bt == TLP_LO || bt == TLP_UL
        S1 = MOI.GreaterThan{Tv}
        throw(MOI.LowerBoundAlreadySet{S1, typeof(s)}(vidx))
    elseif bt == TLP_FX
        S1 = MOI.EqualTo{Tv}
        throw(MOI.LowerBoundAlreadySet{S1, typeof(s)}(vidx))
    elseif bt == TLP_RG
        S1 = MOI.Interval{Tv}
        throw(MOI.LowerBoundAlreadySet{S1, typeof(s)}(vidx))
    end

    # Upper-bound can be added now
    # Two cases: TLP_UP => TLP_UL; TLP_FR => TLP_UP
    if bt == TLP_UP
        # update bounds
        var.dat.lb = s.lower
        # update bound key
        var.dat.bt = TLP_UL
    elseif bt == TLP_FR
        var.dat.lb = s.lower
        var.dat.bt = TLP_LO
    else
        error("Unsupported bound type: $bt")
    end
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}(vidx.uuid)
    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.EqualTo{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    MOI.is_valid(m, f.variable) || throw(MOI.InvalidIndex(f.variable))
    
    # identify which variable
    vidx = VarId(f.variable.value)
    var = m.inner.pbdata_raw.vars[vidx]

    # Get current bounds
    bt, lb, ub = get_bounds(var)

    # Check if upper bound already exists
    if bt == TLP_LO || bt == TLP_UL
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_UP
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_FX
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_RG
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Tv}, typeof(s)}(vidx))
    end

    # Upper-bound can be added now
    # Two cases: TLP_UP => TLP_UL; TLP_FR => TLP_UP
    set_bounds!(var, s.value, s.value)
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}(vidx.uuid)
    return cidx
end

function MOI.add_constraint(
    m::Optimizer{Tv},
    f::MOI.SingleVariable,
    s::MOI.Interval{T}
) where{T<:Real, Tv<:Real}

    # Check that variable exists
    MOI.is_valid(m, f.variable) || throw(MOI.InvalidIndex(f.variable))
    
    # identify which variable
    vidx = VarId(f.variable.value)
    var = m.inner.pbdata_raw.vars[vidx]

    # Get current bounds
    bt, lb, ub = get_bounds(var)

    # Check if upper bound already exists
    if bt == TLP_LO || bt == TLP_UL
        throw(MOI.LowerBoundAlreadySet{MOI.GreaterThan{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_UP
        throw(MOI.UpperBoundAlreadySet{MOI.LessThan{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_FX
        throw(MOI.LowerBoundAlreadySet{MOI.EqualTo{Tv}, typeof(s)}(vidx))
    elseif bt == TLP_RG
        throw(MOI.LowerBoundAlreadySet{MOI.Interval{Tv}, typeof(s)}(vidx))
    end

    # Set bounds on variable
    set_bounds!(var, s.lower, s.upper)
    cidx = MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}(vidx.uuid)
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

    # Create new constraint
    colids = [
        VarId(t.variable_index.value)
        for t in f.terms
    ]
    colvals = [
        t.coefficient
        for t in f.terms
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
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}
) where{Tv<:Real}
    # get index of variable
    vidx = VarId(c.value)

    # Check if variable exists in model
    if haskey(m.inner.pbdata_raw.vars, vidx)
        bt, lb, ub = get_var_bounds(m.inner, vidx)
        return bt == TLP_UP || bt == TLP_UL
    end

    return false
end

function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}
) where{Tv<:Real}
    # get index of variable
    vidx = VarId(c.value)

    # Check if variable exists in model
    if haskey(m.inner.pbdata_raw.vars, vidx)
        bt, lb, ub = get_var_bounds(m.inner, vidx)
        return bt == TLP_LO || bt == TLP_UL
    end

    return false
end
    
function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Tv}}
) where{Tv<:Real}
    # get index of variable
    vidx = VarId(c.value)

    # Check if variable exists in model
    if haskey(m.inner.pbdata_raw.vars, vidx)
        bt, lb, ub = get_var_bounds(m.inner, vidx)
        return bt == TLP_RG
    end

    return false
end

function MOI.is_valid(
    m::Optimizer{Tv},
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.EqualTo{Tv}}
) where{Tv<:Real}
    # get index of variable
    vidx = VarId(c.value)

    # Check if variable exists in model
    if haskey(m.inner.pbdata_raw.vars, vidx)
        bt, lb, ub = get_var_bounds(m.inner, vidx)
        return bt == TLP_FX
    end

    return false
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

        if (bt == TLP_UP || bt == TLP_UL) && S == MOI.LessThan{Tv}
            return MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}(cidx.uuid)
        elseif (bt == TLP_LO || bt == TLP_UL) && S == MOI.GreaterThan{Tv}
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
    m::Optimizer{Tv}, ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)

    return MOI.get(m, MOI.VariablePrimal(), MOI.VariableIndex(c.value))
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintPrimal,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    return MOI.Utilities.get_fallback(m, MOI.ConstraintPrimal(), c)
end

# =============================================
#   ConstraintDual
# =============================================
function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.LessThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable index
    vidx = VarId(c.value)

    # Get variable bounds
    bt, lb, ub = get_var_bounds(m.inner, vidx)
    # Get index in standard form
    vind = m.inner.pbdata_std.var2idx[vidx]

    if bt == TLP_UP
        # x_i <= u was transformed into -x_i >= -u
        # Therefore, we need to return -s_i

        # Get corresponding dual
        s = m.inner.solver.pt.s[vind] / m.inner.solver.pt.t
        return -s
    elseif bt == TLP_UL
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
    m::Optimizer{Tv}, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, MOI.GreaterThan{Tv}}
) where{Tv<:Real}
    MOI.throw_if_not_valid(m, c)  # Sanity check

    # Get variable index
    vidx = VarId(c.value)
    
    # Get index in standard form
    vind = m.inner.pbdata_std.var2idx[vidx]

    # Get corresponding dual
    s = m.inner.solver.pt.s[vind] / m.inner.solver.pt.t
    return s
end

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:Union{MOI.EqualTo{Tv}, MOI.Interval{Tv}}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

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
    m::Optimizer{Tv}, ::MOI.ConstraintDual,
    c::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Tv}, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    MOI.throw_if_not_valid(m, c)  # Sanity check

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


# ==============================================================================
#           V. Objective
# ==============================================================================
    
    # =============================================
    #   V.1 Supported objectives
    # =============================================  
function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}
) where {Tv<:Real}
    return true
end

    # =============================================
    #   V.2 Set objective function
    # =============================================  
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
