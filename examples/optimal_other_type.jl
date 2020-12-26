using Tulip
import MathOptInterface
const MOI = MathOptInterface
using Test

const T = Float32

lp = Tulip.Optimizer{T}()

# Create variables
x = MOI.add_variable(lp)
y = MOI.add_variable(lp)

# Set variable bounds
MOI.add_constraint(lp, MOI.SingleVariable(x), MOI.GreaterThan(T(0)))  # x >= 0
MOI.add_constraint(lp, MOI.SingleVariable(y), MOI.GreaterThan(T(0)))  # y >= 0

# Add constraints
row1 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(T[1.0, -1.0], [x, y]), T(0)),
    MOI.GreaterThan(T(-2))
)
row2 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(T[2.0, -1.0], [x, y]), T(0)),
    MOI.LessThan(T(4))
)
row3 = MOI.add_constraint(lp,
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(T[1.0,  2.0], [x, y]), T(0)),
    MOI.LessThan(T(7))
) 

# Set the objective
MOI.set(lp,
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float32}}(),
    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(T[-2.0, -1.0], [x, y]), T(0))
)
MOI.set(lp, MOI.ObjectiveSense(), MOI.MIN_SENSE)

MOI.optimize!(lp)

objval = MOI.get(lp, MOI.ObjectiveValue())
x_ = MOI.get(lp, MOI.VariablePrimal(), x)
y_ = MOI.get(lp, MOI.VariablePrimal(), y)

@test objval ≈ -8
@test x_ ≈ 3
@test y_ ≈ 2
@test objval isa Float32
@test x_ isa Float32
@test y_ isa Float32
