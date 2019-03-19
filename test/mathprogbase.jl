import MathProgBase

const MPB = MathProgBase

solver_tlp = Tulip.TulipSolver(verbose=1, time_limit=10, algo=1)

m, n, c, b, A, collb, colub, ranges = Tulip.readmps(testdir*"/../dat/netlib/AFIRO.SIF");

# instanciate model
model_tlp = MPB.LinearQuadraticModel(solver_tlp)
@test typeof(model_tlp) <: MPB.AbstractLinearQuadraticModel

MPB.loadproblem!(model_tlp, A, collb, colub, c, b, b, :Min)

@test MPB.getvarLB(model_tlp) == collb
@test MPB.getvarUB(model_tlp) == colub
@test MPB.getobj(model_tlp) == c
@test MPB.getconstrmatrix(model_tlp) == A


MPB.optimize!(model_tlp)


# get solution


# # TODO: Run MathProgBase tests instead
# solver_tlp = Tulip.TulipSolver(verbose=0, time_limit=100.0)
# include(joinpath(dirname(pathof(MathProgBase))[1:(end-3)],"test","linprog.jl"))
# linprogtest(solver_tlp)
# MathProgBase.setparameters!(solver_tlp, Silent=true, TimeLimit=100.0)
# include(joinpath(dirname(pathof(MathProgBase))[1:(end-3)],"test","linproginterface.jl"))
# linprogsolvertest(solver_tlp)