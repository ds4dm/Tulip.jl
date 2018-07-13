include("../src/Tulip.jl")

import MathProgBase
const MPB = MathProgBase

import CPLEX: CplexSolver

function generate_instance(m, n, R; seed=0)

    srand(seed)
    u = rand(m, n*R)

    A = Tulip.LinearAlgebra.DenseBlockAngular(
        [u[:, (1+(r-1)*n):(r*n)] for r in 1:R],
        hcat(eye(m), -eye(m))
    )
    A_ = vcat(
        hcat(kron(speye(R), ones(1, n)), spzeros(R, 2*m)),
        hcat(u, speye(m), -speye(m))
    )

    b = vcat(ones(R), rand(m))
    c = vcat(rand(n*R), 10^4*ones(2*m))

    return A, A_, b, c
end

A, A_, b, c = generate_instance(2, 2, 2)
model = Tulip.Model(A, b, c)
model.env.output_level = 0
Tulip.optimize!(model)

function run_benchmark(m, n, R)
    A, A_, b, c = generate_instance(m, n, R)
    collb = zeros(n*R+2*m)
    colub = Inf * ones(n*R+2*m)
    
    # solve with CPLEX
    solver_cpx = CplexSolver(
        CPX_PARAM_THREADS=1,  # single thread
        CPX_PARAM_LPMETHOD=4, # use barrier algorithm
        CPX_PARAM_PREIND=0,  # no presolve
        CPX_PARAM_SOLUTIONTYPE=2  # no crossover
    )
    model_cpx = MPB.LinearQuadraticModel(solver_cpx)
    MPB.loadproblem!(model_cpx, A_, collb, colub, c, b, b, :Min)
    @time MPB.optimize!(model_cpx)
    
    # solve with Tulip
    solver_tlp = Tulip.TulipSolver(output_level=1, time_limit=120)
    model_tlp = MPB.LinearQuadraticModel(solver_tlp)
    MPB.loadproblem!(model_tlp, A, collb, colub, c, b, b, :Min)
    @time MPB.optimize!(model_tlp)

    println(MPB.getobjval(model_cpx))
    println(MPB.getobjval(model_tlp))
    return nothing
end