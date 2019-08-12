import Base:
    copy,
    getindex, setindex!

include("params.jl")

mutable struct TulipEnv
    #=======================================================
        Algorithmic features
    =======================================================#

    algo::IntParam   # Which interior-point algorithm
    matrix_type::Type  # Type of constraint matrix


    #=======================================================
        Stopping criteria & tolerances
    =======================================================#

    barrier_iter_max::IntParam          # Maximum number of barrier iterations
    time_limit::FloatParam              # Time limit (in seconds)

    barrier_tol_pfeas::FloatParam       # Primal feasibility tolerance
    barrier_tol_dfeas::FloatParam       # Dual feasibility tolerance
    barrier_tol_conv::FloatParam        # Optimality gap tolerance
    barrier_tol_infeas::FloatParam      # Infeasibility tolerance

    beta1::FloatParam
    beta2::FloatParam
    beta3::FloatParam
    beta4::FloatParam

    barrier_max_num_cor::IntParam       # Max number of centrality corrections


    #=======================================================
        Other parameters
    =======================================================#

    verbose::IntParam           # 0 means no output, 1 means normal
    
    # create environment with default values
    # user can over-ride these values afterwards
    function TulipEnv()
        env = new()

        env.algo = RealParam(:algo, 1, 0, 1)
        env.matrix_type = SparseMatrixCSC
        
        env.verbose = RealParam(:verbose, 0, 0, 1)
        env.barrier_iter_max = RealParam(:barrier_iter_max, 100, 0, typemax(Int64))
        env.time_limit = RealParam(:time_limit, Inf, 0.0, Inf)

        env.barrier_tol_pfeas  = RealParam(:barrier_tol_pfeas,  1e-8, 0.0, 1.0)
        env.barrier_tol_dfeas  = RealParam(:barrier_tol_dfeas,  1e-8, 0.0, 1.0)
        env.barrier_tol_conv   = RealParam(:barrier_tol_conv,   1e-8, 0.0, 1.0)
        env.barrier_tol_infeas = RealParam(:barrier_tol_infeas, 1e-8, 0.0, 1.0)

        env.beta1 = RealParam(:beta1, 1e-1,   0.0, 1.0)
        env.beta2 = RealParam(:beta2, 1e-8,   0.0, 1.0)
        env.beta3 = RealParam(:beta3, 0.9999, 0.0, 1.0)
        env.beta4 = RealParam(:beta4, 1e-1,   0.0, 1.0)

        env.barrier_max_num_cor = RealParam(:barrier_max_num_cor, 5, 0, typemax(Int64))

        return env
    end

end

"""
    copy(env)

Copy environmenment.
"""
function copy(env::TulipEnv)
    env_ = TulipEnv()
    for s in fieldnames(TulipEnv)
        p = copy(Core.getfield(env, s))
        Core.setfield!(env_, s, p)
    end
    return env_
end

"""
    getindex(env::TulipEnv, p::Symbol)

Retrieve parameter value. Raises an error if parameter does not exist.

    getindex(env::TulipEnv, p::String)

"""
getindex(env::TulipEnv, p::Symbol) = get_param_value(Core.getfield(env, p))
getindex(env::TulipEnv, param::String) = getindex(env, Symbol(param))

# import Base.getproperty
# function getproperty(env::TulipEnv, name::Symbol)
#     p = Core.getfield(env, name)
#     if isa(p, AbstractParam)
#         return p.val
#     else
#         return p
#     end
# end

"""
    setindex!(env, v, p)

Set value of given parameter to v. Raises an error if parameter does not exist,
or if incorrect value.
"""
setindex!(env::TulipEnv, v, p::Symbol) = set_param_value!(Core.getfield(env, p), v)
setindex!(env::TulipEnv, v, p::String) = set_param_value!(Core.getfield(env, Symbol(p)), v)

"""
    set_param!(env; kwargs...)

Set multiple parameters at the same time.
"""
function set_param!(env; kwargs...)
    for (p, v) in kwargs
        set_param_value!(Core.getfield(env, p), v)
    end
end

"""
    reset!(env)

Reset all parameters' values to default
"""
function reset!(env::TulipEnv)
    for p in fieldnames(TulipEnv)
        set_param_default!(Core.getfield(env, p))
    end
    return nothing
end