print("\ttst_basic.jl")

# generate random instance
srand(0)
m = 2^4
n = 2^4
A = hcat(sprand(m, n, 0.1), speye(m))
u = 10.0 * sprand(n, 1.0)
b = 1.0 + rand(m)
c = vcat(1.0 - 2.0 * rand(n), zeros(m))

model = Tulip.Model(A, b, c, u.nzind, u.nzval)
status = Tulip.solve!(model, verbose=0, tol=10.0^-8)

@test status == :Optimal

objval = dot(model.sol.x, model.c)
const OBJVAL = -15.989485968940144
@test abs(objval - OBJVAL) <= 10.0^-8

println("\tPassed.")