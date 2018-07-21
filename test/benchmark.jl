using BenchmarkTools
Base.BLAS.set_num_threads(1)

import MathProgBase
const MPB = MathProgBase

import Gurobi: GurobiSolver
import Clp: ClpSolver
import GLPKMathProgInterface: GLPKSolverLP
import CPLEX: CplexSolver
import Mosek: MosekSolver
import OSQP

import Tulip


function generate_instance(m, n, R; seed=0)
    
    srand(seed)

    u = [1.0 - 2.0 * rand(m, n) for _ in 1:R]
    for r in 1:R
        u[r][:, 1] = 0.0
    end
    colptr = [1+n*(r-1) for r in 1:(R+1)]
    A = Tulip.LinearAlgebra.DenseBlockAngular(u)
    b = vcat(ones(R), zeros(m))
    c = rand(n*R)
    colub_ind = Vector{Int}(0,)
    colub_val = Vector{Float64}(0,)

    senses = ['=' for i=1:(m+R)]
    u_ = Inf * ones(n*R)

    e_ = kron(speye(R), sparse(ones(1, n)));
    A_ = vcat(e_, hcat(u...));
    
    return (A, A_, b, c, colub_ind, colub_val, senses)
end

function run_benchmark(A, A_, c, rowlb, rowub, collb, colub, sense, solver::MPB.AbstractMathProgSolver, name::Symbol)
    # println(name)
    
    # load model
    model = MPB.LinearQuadraticModel(solver)
    if name == :TULIP || name == :TLP
        MPB.loadproblem!(model, A, collb, colub, c, rowlb, rowlb, sense)
    else
        MPB.loadproblem!(model, A_, collb, colub, c, rowlb, rowlb, sense)
    end
    
    t = @elapsed MPB.optimize!(model)
    
    # get results
    s = MPB.status(model)
    obj = MPB.getobjval(model)
    
    return (name, s, obj, t)
end

"""
    generate_solver()

"""
generate_solver(s::Symbol;time_limit=300, nthreads::Int=1, verbose::Int=0) = generate_solver(Val{s}; time_limit=time_limit, nthreads=nthreads, verbose=verbose)
function generate_solver(::Type{Val{:CPX}};time_limit=300, nthreads::Int=1, verbose::Int=0)
    CplexSolver(
        CPX_PARAM_SCRIND=verbose,
        CPX_PARAM_LPMETHOD=4,
        CPX_PARAM_PREIND=0,
        CPX_PARAM_THREADS=nthreads,
        CPX_PARAM_TILIM=time_limit
    )
end

function generate_solver(::Type{Val{:CLP}};time_limit=300, nthreads::Int=1, verbose::Int=0)
    ClpSolver(
        LogLevel=verbose,
        SolveType=4,
        PresolveType=-1,
        MaximumSeconds=time_limit
    )
end

function generate_solver(::Type{Val{:GLPK}};time_limit=300, nthreads::Int=1, verbose::Int=0)
    GLPKSolverLP(presolve=false, method=:InteriorPoint)
end

function generate_solver(::Type{Val{:GRB}};time_limit=300, nthreads::Int=1, verbose::Int=0)
    GurobiSolver(
        OutputFlag=verbose,
        Method=2,
        Presolve=0,
        Threads=nthreads,
        Crossover=0,
        TimeLimit=time_limit
    )
end

function generate_solver(::Type{Val{:OSQP}};time_limit=300, nthreads=1, verbose::Int=0)
    OSQP.OSQPMathProgBaseInterface.OSQPSolver(
        verbose=verbose,
        time_limit=time_limit
    )
end

function generate_solver(::Type{Val{:MSK}};time_limit=300, nthreads=1, verbose::Int=0)
    MosekSolver(
        MSK_IPAR_OPTIMIZER=4,
        MSK_IPAR_PRESOLVE_USE=0,
        MSK_IPAR_NUM_THREADS=nthreads,
        MSK_IPAR_INTPNT_BASIS=0,
        MSK_DPAR_OPTIMIZER_MAX_TIME=time_limit,
        MSK_IPAR_LOG=verbose
    )
end

function generate_solver(::Type{Val{:TLP}};time_limit=300, nthreads=1, verbose::Int=0)
    Tulip.TulipSolver(
        algo=1,
        output_level=verbose,
        time_limit=time_limit
    )
end

generate_solvers(solver_list;time_limit=300, nthreads=1, verbose::Int=0) = [
    (generate_solver(name, time_limit=time_limit, nthreads=nthreads, verbose=verbose), name)
    for name in solver_list
]

#====================================
    Benchmarks
=====================================#

# full benchmark
function benchmark(
    size_instances::Vector{Tuple{Int64, Int64, Int64}},
    solvers::Vector{Tuple{MPB.AbstractMathProgSolver, Symbol}},
    seeds::Vector{Int64};
    resultsfile::String=""
)

    results = Dict{Tuple{Int64,Int64,Int64,Int64,Symbol}, Tuple{Float64,Float64}}()
    # write results to file
    # res_filename = "benchmark_"*Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")*".csv"
    if resultsfile != ""
        open(resultsfile, "w") do f
            write(f, "m, n, R, s, solver, time\n")
        end
    end

    for (m, n, R) in size_instances
        # run 5 seeds
        println((m, n, R))
        for seed in seeds
            # generate instance
            (A, A_, b, c, colub_ind, colub_val, senses) = generate_instance(m, n, R, seed=seed)
            for (solver, name) in solvers
                # println(name)
                res = run_benchmark(A, A_, c, b, b, zeros(n*R), Inf*ones(n*R), :Min, solver, name)
                results[(m, n, R, seed, name)] = (res[4], res[3])
                if resultsfile != ""
                    open(resultsfile, "a") do f
                        write(f, "$(m), $(n), $(R), $(seed), $(name), $(res[4])\n")
                    end
                end
            end
        end
    end

    return results
end