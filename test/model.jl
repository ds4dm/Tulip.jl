#= Create and solve the following model:
    min  x1 + 2*x2
    s.t. x1 + x2 = 2
            x2 < 1
            x1, x2 >0

The solution to this problem is:
    x1 = 2.0, x2 = 0.0
    y = 1.0
    s1 = 0.0, s2 = 1.0
    w2 = 1.0
    z2 = 0.0
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
@test abs(2.0 - Tulip.getobjectivevalue(model)) <= 10.0^-8
@test abs(2.0 - Tulip.getdualbound(model)) <= 10.0^-8
@test Tulip.getobjectivedualgap(model) <= 10.0^-8
@test Tulip.getnumbarrieriter(model) <= 100

x_ = Tulip.getsolution(model)
y_ = Tulip.getconstrduals(model)
s_ = Tulip.getreducedcosts(model)

@test norm(x_ - [2.0, 0.0], Inf) <= 10.0^-8
@test norm(y_ - [1.0], Inf) <= 10.0^-8
@test norm(s_ - [0.0, 1.0], Inf) <= 10.0^-8