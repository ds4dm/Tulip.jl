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
    elseif st == Sln_Optimal || Sln_FeasiblePoint
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
    var_bounds_name::Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}

    function Optimizer{Tv}() where{Tv<:Real}
        return new{Tv}(
            Model{Tv}(),
            Dict{MOI.ConstraintIndex{MOI.SingleVariable, <:SCALAR_SETS{Tv}}, String}()
        )
    end
end

# TODO: set parameters via kw arguments
Optimizer() = Optimizer{Float64}()

MOI.empty!(m::Optimizer) = empty!(m.inner)

MOI.is_empty(m::Optimizer) = is_empty(m.inner)

MOI.optimize!(m::Optimizer) = optimize!(m.inner)

# ==============================================================================
#           I. Optimizer attributes
# ==============================================================================

SUPPORTED_OPTIMIZER_ATTR = Union{
    MOI.SolverName,
    MOI.Silent,
    MOI.TimeLimitSec
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_OPTIMIZER_ATTR} = true

MOI.get(::Optimizer, ::MOI.SolverName) = "Tulip"

MOI.get(m::Optimizer, ::MOI.Silent) = iszero(m.inner.env.verbose)

function MOI.set(m::Optimizer, ::MOI.Silent, flag::Bool)
    m.inner.env.verbose = 1 - flag
    return nothing
end

MOI.get(m::Optimizer, ::MOI.TimeLimitSec) = m.inner.env.time_limit

function MOI.set(m::Optimizer, ::MOI.TimeLimitSec, t)
    m.inner.env.time_limit = t
    return nothing
end

# ==============================================================================
#           II. Model attributes
# ==============================================================================

SUPPORTED_MODEL_ATTR = Union{
    MOI.Name,
    MOI.ObjectiveSense,
    MOI.NumberOfVariables,
    MOI.ListOfVariableIndices,
    MOI.ListOfConstraintIndices,
    MOI.NumberOfConstraints,
    # ListOfConstraints,  # TODO
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction},
    MOI.ObjectiveFunctionType,
    MOI.ObjectiveValue,
    MOI.DualObjectiveValue,
    # MOI.ObjectiveBound,  # This is generally used for MIP models
    MOI.RelativeGap,
    # MOI.SolveTime,  # TODO
    MOI.SimplexIterations,
    MOI.BarrierIterations,
    MOI.RawSolver,
    MOI.ResultCount,
    MOI.TerminationStatus,
    MOI.PrimalStatus,
    MOI.DualStatus
}

MOI.supports(::Optimizer, ::A) where{A<:SUPPORTED_MODEL_ATTR} = true

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

MOI.get(m::Optimizer, ::MOI.Name) = m.inner.name

MOI.get(m::Optimizer, ::MOI.NumberOfVariables) = get_num_var(m.inner)

function MOI.get(m::Optimizer, ::MOI.ListOfVariableIndices)
    return [
        MOI.VariableIndex(vidx.uuid)
        for vidx in keys(m.pbdata_raw.vars)
    ]
end

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

MOI.get(
    m::Optimizer{Tv}, ::MOI.ObjectiveFunctionType
) where{Tv<:Real} = MOI.ScalarAffineFunction{Tv}

function MOI.get(m::Optimizer{Tv}, ::MOI.ObjectiveValue) where{Tv<:Real}
    # TODO: dispatch a function call on m.inner
    return m.inner.solver.primal_bound_scaled
end

function MOI.get(m::Optimizer{Tv}, ::MOI.DualObjectiveValue) where{Tv<:Real}
    # TODO: dispatch a function call on m.inner
    return m.inner.solver.dual_bound_scaled
end

function MOI.get(m::Optimizer{Tv}, ::MOI.RelativeGap) where{Tv<:Real}
    # TODO: dispatch a function call on m.inner
    zp = m.inner.solver.primal_bound_scaled
    zd = m.inner.solver.dual_bound_scaled
    return (abs(zp - zd) / (Tv(1e-10) + abs(zd)))
end

# TODO: MOI.get(m::Optimizer, ::MOI.SolveTime)

MOI.get(::Optimizer, ::MOI.SimplexIterations) = 0

# TODO: dispatch function call on m.inner
MOI.get(m::Optimizer, ::MOI.BarrierIterations) = m.inner.solver.niter

MOI.get(m::Optimizer, ::MOI.RawSolver) = m.inner

MOI.get(::Optimizer, ::MOI.ResultCount) = 1

function MOI.get(m::Optimizer, ::MOI.TerminationStatus)
    # TODO: dispatch function call on m.inner
    isa(m.inner.solver, Nothing) && return MOI.OPTIMIZE_NOT_CALLED

    st = m.inner.solver.solver_status
    
    return MOITerminationStatus(st)
end

MOI.get(m::Optimizer, ::MOI.PrimalStatus) = MOISolutionStatus(m.inner.solver.primal_status)

MOI.get(m::Optimizer, ::MOI.DualStatus) = MOISolutionStatus(m.inner.solver.dual_status)


# ==============================================================================
#           III. Variables
# ==============================================================================

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

# TODO: delete variable (and multiple variables)
# function MOI.delete(m::Optimizer{Tv}, v::MOI.VariableIndex) where{Tv<:Real} end

# TODO: query variable by name
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

function MOI.get(
    m::Optimizer{Tv}, ::MOI.VariablePrimal, x::MOI.VariableIndex
) where{Tv<:Real}
    MOI.is_valid(m, x) || throw(MOI.InvalidIndex(x))
    return get_value(m.inner, VarId(x.value))
end


# ==============================================================================
#           IV. Constraints
# ==============================================================================

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
    # MOI.ConstraintFunction,  # TODO
    MOI.ConstraintSet
}

MOI.supports(::Optimizer, ::A, ::Type{<:MOI.ConstraintIndex}) where{A<:SUPPORTED_CONSTR_ATTR} = true


# Variable bounds
# TODO: change MOI convention to handle explicit bounds
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.SingleVariable}, ::Type{S}
) where {Tv<:Real, S<:SCALAR_SETS{Tv}}
    return true
end

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

# Linear constraints
function MOI.supports_constraint(
    ::Optimizer, ::Type{MOI.ScalarAffineFunction{Tv}}, ::Type{S}
) where {Tv<:Real, S<:SCALAR_SETS{Tv}}
    return true
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

function MOI.get(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S}
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}
    
    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Get constraint name for dict directly
    if !haskey(m.var_bounds_name, c)
        return ""
    else
        return m.var_bounds_name[c]
    end
end

function MOI.set(
    m::Optimizer{Tv}, ::MOI.ConstraintName,
    c::MOI.ConstraintIndex{MOI.SingleVariable, S},
    name::String
) where{Tv<:Real, S<:SCALAR_SETS{Tv}}

    MOI.is_valid(m, c) || throw(MOI.InvalidIndex(c))

    # Set constraint name in dict directly
    m.var_bounds_name[c] = name
    return nothing
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
    return nothing
end


# ==============================================================================
#           V. Objective
# ==============================================================================

function MOI.supports(
    ::Optimizer,
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Tv}}
) where {Tv<:Real}
    return true
end

function MOI.set(
    m::Optimizer{Tv},
    ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
    f::MOI.ScalarAffineFunction{T}
) where{T<:Real, Tv<:Real}

    # Check that constant term is zero
    iszero(f.constant) || error("Objective constant must be zero.")

    # First, set objective to zero
    for (vidx, var) in m.inner.pbdata_raw.vars
        set_obj_coeff!(var, zero(Tv))
    end

    # Update new coefficients
    for t in f.terms
        var = m.inner.pbdata_raw.vars[VarId(t.variable_index.value)]
        set_obj_coeff!(var, t.coefficient)
    end

    return nothing
end