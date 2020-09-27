# Toy example

Tulip can be accessed in 3 ways:
through [JuMP](https://github.com/jump-dev/JuMP.jl),
through [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl),
or directly.

This tutorial illustrates, for each case, how to build a model, solve it,
and query the solution value.
In all cases, we consider the small LP
```math
\begin{array}{rrrl}
    (LP) \ \ \ 
    \displaystyle Z^{*} = \min_{x, y} & -2x & - y\\
    s.t.
    &  x & - y & \geq -2\\
    & 2x &-  y & \leq  4\\
    &  x &+ 2y & \leq  7\\
    &  x,&   y & \geq  0\\
\end{array}
```
whose optimal value and solution are ``Z^{*} = -8`` and ``(x^{*}, y^{*}) = (3, 2)``.

## JuMP

```jldoctest; output = false
using Printf
using JuMP
import Tulip

# Instantiate JuMP model
lp = Model(Tulip.Optimizer)

# Create variables
@variable(lp, x >= 0)
@variable(lp, y >= 0)

# Add constraints
@constraint(lp, row1, x - y >= -2)
@constraint(lp, row2, 2*x - y <= 4)
@constraint(lp, row3, x + 2*y <= 7)

# Set the objective
@objective(lp, Min, -2*x - y)

# Set some parameters
set_optimizer_attribute(lp, "OutputLevel", 0)  # disable output
set_optimizer_attribute(lp, "Presolve", 0)     # disable presolve

# Solve the problem
optimize!(lp)

# Check termination status
st = termination_status(lp)
println("Termination status: $st")

# Query solution value
objval = objective_value(lp)
x_ = value(x)
y_ = value(y)

@printf "Z* = %.4f\n" objval
@printf "x* = %.4f\n" x_
@printf "y* = %.4f\n" y_

# output

Termination status: OPTIMAL
Z* = -8.0000
x* = 3.0000
y* = 2.0000
```

## MOI

```jldoctest; output = false
using Printf

import MathOptInterface
const MOI = MathOptInterface

import Tulip

lp = Tulip.Optimizer{Float64}()

# Create variables
x = MOI.add_variable(lp)
y = MOI.add_variable(lp)

# Set variable bounds
MOI.add_constraint(lp, MOI.SingleVariable(x), MOI.GreaterThan(0.0))  # x >= 0
MOI.add_constraint(lp, MOI.SingleVariable(y), MOI.GreaterThan(0.0))  # y >= 0

# Add constraints
row1 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, -1.0], [x, y]), 0.0),
    MOI.GreaterThan(-2.0)
)
row2 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0, -1.0], [x, y]), 0.0),
    MOI.LessThan(4.0)
)
row3 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0,  2.0], [x, y]), 0.0),
    MOI.LessThan(7.0)
) 

# Set the objective
MOI.set(lp,
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-2.0, -1.0], [x, y]), 0.0)
)
MOI.set(lp, MOI.ObjectiveSense(), MOI.MIN_SENSE)

# Set some parameters
MOI.set(lp, MOI.Silent(), true)               # disable output
MOI.set(lp, MOI.RawParameter("Presolve"), 0)  # disable presolve

# Solve the problem
MOI.optimize!(lp)

# Check status
st = MOI.get(lp, MOI.TerminationStatus())
println("Termination status: $st")

# Query solution value
objval = MOI.get(lp, MOI.ObjectiveValue())
x_ = MOI.get(lp, MOI.VariablePrimal(), x)
y_ = MOI.get(lp, MOI.VariablePrimal(), y)

@printf "Z* = %.4f\n" objval
@printf "x* = %.4f\n" x_
@printf "y* = %.4f\n" y_

# output

Termination status: OPTIMAL
Z* = -8.0000
x* = 3.0000
y* = 2.0000
```

## Tulip

!!! warning
    Tulip's low-level API should not be considered stable nor complete.
    The recommended way to use Tulip is through JuMP/MOI as shown above.


```jldoctest; output = false
using Printf
import Tulip

# Instantiate Tulip object
lp = Tulip.Model{Float64}()
pb = lp.pbdata  # Internal problem data

# Create variables
x = Tulip.add_variable!(pb, Int[], Float64[], -2.0, 0.0, Inf, "x")
y = Tulip.add_variable!(pb, Int[], Float64[], -1.0, 0.0, Inf, "y")

# Add constraints
row1 = Tulip.add_constraint!(pb, [x, y], [1.0, -1.0], -2.0, Inf, "row1")
row2 = Tulip.add_constraint!(pb, [x, y], [2.0, -1.0], -Inf, 4.0, "row2")
row3 = Tulip.add_constraint!(pb, [x, y], [1.0,  2.0], -Inf, 7.0, "row3")

# Set the objective
# Nothing to do here as objective is already declared

# Set some parameters
Tulip.set_parameter(lp, "OutputLevel", 0)  # disable output
Tulip.set_parameter(lp, "Presolve", 0)     # disable presolve

# Solve the problem
Tulip.optimize!(lp)

# Check termination status
st = Tulip.get_attribute(lp, Tulip.Status())
println("Termination status: $st")

# Query solution value
objval = Tulip.get_attribute(lp, Tulip.ObjectiveValue())
x_ = lp.solution.x[x]
y_ = lp.solution.x[y]

@printf "Z* = %.4f\n" objval
@printf "x* = %.4f\n" x_
@printf "y* = %.4f\n" y_

# output

Termination status: Trm_Optimal
Z* = -8.0000
x* = 3.0000
y* = 2.0000
```