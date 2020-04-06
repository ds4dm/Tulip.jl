```@meta
CurrentModule = Tulip
```

# Attribute reference

Attributes are queried using [`get_attribute`](@ref) and set using [`set_attribute`](@ref).

## Model attributes

| Name                              | Type      | Description                   
|:----------------------------------|:----------|:------------------------------
| [`ModelName`](@ref)               | `String`  | Name of the model
| [`NumberOfConstraints`](@ref)     | `Int`     | Number of constraints in the model
| [`NumberOfVariables`](@ref)       | `Int`     | Number of variables in the model
| [`ObjectiveValue`](@ref)          | `Tv`      | Objective value of the current primal solution
| [`DualObjectiveValue`](@ref)      | `Tv`      | Objective value of the current dual solution
| [`ObjectiveConstant`](@ref)       | `Tv`      | Value of the objective constant
| [`ObjectiveSense`](@ref)          |           | Optimization sense
| [`Status`](@ref)                  |           | Model status
| [`BarrierIterations`](@ref)       | `Int`     | Number of barrier iterations
| [`SolutionTime`](@ref)            | `Float64` | Solution time, in seconds

## Variable attributes

| Name                               | Type     | Description                   
|:-----------------------------------|:---------|:------------------------------
| [`VariableLowerBound`](@ref)       | `Tv`     | Variable lower bound 
| [`VariableUpperBound`](@ref)       | `Tv`     | Variable upper bound 
| [`VariableObjectiveCoeff`](@ref)   | `Tv`     | Variable objective coefficient 
| [`VariableName`](@ref)             | `String` | Variable name 

## Constraint attributes

| Name                               | Type     | Description                   
|:-----------------------------------|:---------|:------------------------------
| [`ConstraintLowerBound`](@ref)     | `Tv`     | Constraint lower bound 
| [`ConstraintUpperBound`](@ref)     | `Tv`     | Constraint upper bound
| [`ConstraintName`](@ref)           | `String` | Constraint name 



## Reference

### Model attributes

```@autodocs
Modules = [Tulip]
Pages = ["src/Core/attributes.jl"]
Filter = t -> typeof(t) === DataType && t <: Tulip.AbstractModelAttribute
```

### Variable attributes

```@autodocs
Modules = [Tulip]
Pages = ["src/Core/attributes.jl"]
Filter = t -> typeof(t) === DataType && t <: Tulip.AbstractVariableAttribute
```

### Constraint attributes

```@autodocs
Modules = [Tulip]
Pages = ["src/Core/attributes.jl"]
Filter = t -> typeof(t) === DataType && t <: Tulip.AbstractConstraintAttribute
```