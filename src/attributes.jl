abstract type AbstractAttribute end

# ==============================================================================
#
#      Model attributes
#
# ==============================================================================
abstract type AbstractModelAttribute <: AbstractAttribute end

"""
    ModelName

The name of the model.

**Type:** `String`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.ModelName(), "lp_example")
Tulip.get_attribute(model, Tulip.ModelName())
```
"""
struct ModelName <: AbstractModelAttribute end

"""
    NumberOfConstraints

Number of constraints in the model.

**Type:** `Int`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.NumberOfConstraints())
```
"""
struct NumberOfConstraints <: AbstractModelAttribute end

"""
    NumberOfVariables

Number of variables in the model.

**Type:** `Int`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.NumberOfVariables())
```
"""
struct NumberOfVariables <: AbstractModelAttribute end

"""
    ObjectiveValue

Objective value of the current primal solution.

**Type:** `T`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.ObjectiveValue())
```
"""
struct ObjectiveValue <: AbstractModelAttribute end

"""
    DualObjectiveValue

Objective value of the current dual solution.

**Type:** `T`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.DualObjectiveValue())
```
"""
struct DualObjectiveValue <: AbstractModelAttribute end

"""
    ObjectiveConstant

Constant objective offset, defaults to zero.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.ObjectiveConstant(), zero(T))
Tulip.get_attribute(model, Tulip.ObjectiveConstant())
```
"""
struct ObjectiveConstant <: AbstractModelAttribute end

"""
    ObjectiveSense

"""
struct ObjectiveSense <: AbstractModelAttribute end

"""
    Status

Model status

### Type:

**Modifiable:** No

### Examples

```julia
Tulip.get(model, Tulip.Status())
```
"""
struct Status <: AbstractModelAttribute end

"""
    BarrierIterations

Number of iterations of the barrier algorithm in the last call.

This number may be zero in the following cases:
* the model has been solved yet
* presolve solved the model
* the initial solution was optimal

**Type:** `Int`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.BarrierIterations())
```
"""
struct BarrierIterations <: AbstractModelAttribute end

"""
    SolutionTime

Total solution time, in seconds.

**Type:** `Float64`

**Modifiable:** No

### Examples

```julia
Tulip.get_attribute(model, Tulip.SolutionTime())
```
"""
struct SolutionTime <: AbstractModelAttribute end


# ==============================================================================
#
#      Variable attributes
#
# ==============================================================================
abstract type AbstractVariableAttribute <: AbstractAttribute end

"""
    VariableLowerBound

Variable lower bound.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.VariableLowerBound(), 1, zero(T))
Tulip.get_attribute(model, Tulip.VariableLowerBound(), 1)
```
"""
struct VariableLowerBound <: AbstractVariableAttribute end

"""
    VariableUpperBound

Variable upper bound.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.VariableUpperBound(), 1, one(T))
Tulip.get_attribute(model, Tulip.VariableUpperBound(), 1)
```
"""
struct VariableUpperBound <: AbstractVariableAttribute end

"""
    VariableObjectiveCoeff

Objective coefficient of the variable.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.VariableObjectiveCoeff(), 1, one(T))
Tulip.get_attribute(model, Tulip.VariableObjectiveCoeff(), 1)
```
"""
struct VariableObjectiveCoeff <: AbstractVariableAttribute end

"""
    VariableName

Name of the variable.

**Type:** `String`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.VariableName(), 1, "x1")
Tulip.get_attribute(model, Tulip.VariableName(), 1)
```
"""
struct VariableName <: AbstractVariableAttribute end


# ==============================================================================
#
#      Constraint attributes
#
# ==============================================================================
abstract type AbstractConstraintAttribute <: AbstractAttribute end

"""
    ConstraintLowerBound

Constraint lower bound.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.ConstraintLowerBound(), 1, zero(T))
Tulip.get_attribute(model, Tulip.ConstraintLowerBound(), 1)
```
"""
struct ConstraintLowerBound <: AbstractConstraintAttribute end

"""
    ConstraintUpperBound

Constraint upper bound.

**Type:** `T`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.ConstraintUpperBound(), 1, one(T))
Tulip.get_attribute(model, Tulip.ConstraintUpperBound(), 1)
```
"""
struct ConstraintUpperBound <: AbstractConstraintAttribute end

"""
    ConstraintName

Name of the constraint.

**Type:** `String`

**Modifiable:** Yes

### Examples

```julia
Tulip.set_attribute(model, Tulip.ConstraintName(), 1, "c1")
Tulip.get_attribute(model, Tulip.ConstraintName(), 1)
```
"""
struct ConstraintName <: AbstractConstraintAttribute end