import MathProgBase

const MPB = MathProgBase

solver_tlp = Tulip.TulipSolver(output_level=0, time_limit=10)

m, n, c, b, A, collb, colub, ranges = Tulip.readmps(testdir*"/../dat/netlib/AFIRO.SIF");

# instanciate model
model_tlp = MPB.LinearQuadraticModel(solver_tlp)
MPB.loadproblem!(model_tlp, A, collb, colub, c, b, b, :Min)