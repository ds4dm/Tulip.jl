print("\tmodel.jl")

#= Create and solve the following model:
    min  x + 2*y
    s.t. x + y = 2
            y < 1
            x, y >0
=#
A = sparse([1.0 1.0])
c = [1.0, 2.0]
b = [2.0]
u = sparse([0.0, 1.0])

m = 1
p = nnz(u)
n = 2

# Create random instance and run checks
model = Tulip.Model(A, b, c, u.nzind, u.nzval)
model.env[:output_level] = 0
Tulip.optimize!(model)

# Low-level interface
@test n == Tulip.getnumvar(model)
@test m == Tulip.getnumconstr(model)
@test [0.0, 0.0] == Tulip.getvarlowerbounds(model)
@test [Inf, 1.0] == Tulip.getvarupperbounds(model)
@test [2.0] == Tulip.getconstrlowerbounds(model)
@test [2.0] == Tulip.getconstrupperbounds(model)
@test [1.0, 2.0] == Tulip.getobjectivecoeffs(model)
@test A == Tulip.getlinearconstrcoeffs(model)
x_ = Tulip.getsolution(model)
y_ = Tulip.getconstrduals(model)
s_ = Tulip.getreducedcosts(model)
@test norm(x_ - [2.0, 0.0], Inf) <= 10.0^-8
@test norm(y_ - [1.0], Inf) <= 10.0^-8
@test norm(s_ - [0.0, 1.0], Inf) <= 10.0^-8

println("\tPassed.")