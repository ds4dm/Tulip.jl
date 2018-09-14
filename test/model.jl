env = Tulip.TulipEnv()
model = Tulip.Model(env)

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

Tulip.addvar!(model, [], 0.0, Inf, 1.0)
Tulip.addvar!(model, [], 0.0, 1.0, 2.0)
Tulip.addconstr!(model, [1.0, 1.0], 2.0, 2.0)

# Tulip.prepross!(model)

A_ = sparse([[1.0 1.0];])

@test model.A == A_
@test model.var_lb == [0.0, 0.0]
@test model.var_ub == [Inf, 1.0]
@test model.constr_lb == [2.0]
@test model.constr_ub == [2.0]

model.env[:verbose] = 1
Tulip.optimize!(model)

# Low-level interface
@test Tulip.getnumvar(model) == 2
@test Tulip.getnumconstr(model) == 1
@test Tulip.getvarlowerbounds(model) == [0.0, 0.0]
@test Tulip.getvarupperbounds(model) == [Inf, 1.0]
@test Tulip.getconstrlowerbounds(model) == [2.0]
@test Tulip.getconstrupperbounds(model) == [2.0]
@test Tulip.getobjectivecoeffs(model) == [1.0, 2.0]
@test Tulip.getlinearconstrcoeffs(model) == A_
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