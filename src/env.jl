import Base.copy

mutable struct Env{Tv<:Real}
    #=======================================================
        Algorithmic features
    =======================================================#

    algo::Int   # Which interior-point algorithm
    matrix_type::Type  # Type of constraint matrix

    ls_backend::TLA.LSBackend
    ls_system::TLA.LinearSystem


    #=======================================================
        Stopping criteria & tolerances
    =======================================================#

    barrier_iter_max::Int          # Maximum number of barrier iterations
    time_limit::Float64            # Time limit (in seconds)

    barrier_tol_pfeas::Tv       # Primal feasibility tolerance
    barrier_tol_dfeas::Tv       # Dual feasibility tolerance
    barrier_tol_conv::Tv        # Optimality gap tolerance
    barrier_tol_infeas::Tv      # Infeasibility tolerance

    beta1::Tv
    beta2::Tv
    beta3::Tv
    beta4::Tv

    barrier_max_num_cor::Int       # Max number of centrality corrections


    #=======================================================
        Other parameters
    =======================================================#

    verbose::Int           # 0 means no output, 1 means normal
    threads::Int           # Number of threads. Should be > 0
    
    # create environment with default values
    # user can over-ride these values afterwards
    function Env{Tv}() where{Tv<:Real}
        env = new()

        env.algo = 1
        env.matrix_type = SparseMatrixCSC

        env.ls_backend = TLPLinearAlgebra.DefaultBackend()
        env.ls_system = TLPLinearAlgebra.DefaultSystem()
        
        env.verbose = 0
        env.threads = 1
        env.barrier_iter_max = 100
        env.time_limit = Inf

        env.barrier_tol_pfeas  = Tv(1e-8)
        env.barrier_tol_dfeas  = Tv(1e-8)
        env.barrier_tol_conv   = Tv(1e-8)
        env.barrier_tol_infeas = Tv(1e-8)

        env.beta1 = Tv(1e-1)
        env.beta2 = Tv(1e-8)
        env.beta3 = Tv(0.9999)
        env.beta4 = Tv(1e-1)

        env.barrier_max_num_cor = 5

        return env
    end

end

"""
    copy(env)

Copy environmenment.
"""
function copy(env::Env{Tv}) where{Tv<:Real}
    env_ = Env{Tv}()
    for s in fieldnames(Env)
        p = Core.getfield(env, s)
        Core.setfield!(env_, s, p)
    end
    return env_
end