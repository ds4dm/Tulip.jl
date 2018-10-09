import Base:
    copy,
    getindex, setindex!

include("params.jl")

mutable struct TulipEnv
    #=======================================================
        Algorithmic features
    =======================================================#

    algo::IntParam   # Which interior-point algorithm


    #=======================================================
        Stopping criteria & tolerances
    =======================================================#

    barrier_iter_max::IntParam       # Maximum number of barrier iterations
    time_limit::FloatParam         # Time limit (in seconds)

    barrier_tol_pfeas::FloatParam   # Primal feasibility tolerance
    barrier_tol_dfeas::FloatParam    # Dual feasibility tolerance
    barrier_tol_conv::FloatParam   # Optimality gap tolerance
    barrier_tol_infeas::FloatParam   # Infeasibility tolerance

    beta1::FloatParam
    beta2::FloatParam
    beta3::FloatParam
    beta4::FloatParam


    #=======================================================
        Other parameters
    =======================================================#

    verbose::IntParam           # 0 means no output, 1 means normal
    
    # create environment with default values
    # user can over-ride these values afterwards
    function TulipEnv()
        env = new()

        env.algo = RealParam(:algo, 1, 0, 1)
        env.verbose = RealParam(:verbose, 0, 0, 1)
        env.barrier_iter_max = RealParam(:barrier_iter_max, 100, 0, typemax(Int64))
        env.time_limit = RealParam(:time_limit, Inf, 0.0, Inf)

        env.barrier_tol_pfeas = RealParam(:barrier_tol_pfeas, 10.0^-8, 0.0, 1.0)
        env.barrier_tol_dfeas = RealParam(:barrier_tol_dfeas, 10.0^-8, 0.0, 1.0)
        env.barrier_tol_conv = RealParam(:barrier_tol_conv, 10.0^-8, 0.0, 1.0)
        env.barrier_tol_infeas = RealParam(:barrier_tol_infeas, 10.0^-8, 0.0, 1.0)

        env.beta1 = RealParam(:beta1, 0.1, 0.0, 1.0)
        env.beta2 = RealParam(:beta2, 10.0^-8, 0.0, 1.0)
        env.beta3 = RealParam(:beta3, 0.9999, 0.0, 1.0)
        env.beta4 = RealParam(:beta4, 0.1, 0.0, 1.0)

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

    getindex(env::TulipEnv, ::Type{Val{:<param>}})
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

getindex(env::TulipEnv, ::Type{Val{:algo}}) = env.algo.val

getindex(env::TulipEnv, ::Type{Val{:verbose}}) = env.verbose.val

getindex(env::TulipEnv, ::Type{Val{:barrier_iter_max}}) = env.barrier_iter_max.val

getindex(env::TulipEnv, ::Type{Val{:time_limit}}) = env.time_limit.val

getindex(env::TulipEnv, ::Type{Val{:barrier_tol_pfeas}}) = env.barrier_tol_pfeas.val

getindex(env::TulipEnv, ::Type{Val{:barrier_tol_dfeas}}) = env.barrier_tol_dfeas.val

getindex(env::TulipEnv, ::Type{Val{:barrier_tol_conv}}) = env.barrier_tol_conv.val

getindex(env::TulipEnv, ::Type{Val{:barrier_tol_infeas}}) = env.barrier_tol_infeas.val

getindex(env::TulipEnv, ::Type{Val{:beta1}}) = env.beta1.val

getindex(env::TulipEnv, ::Type{Val{:beta2}}) = env.beta2.val

getindex(env::TulipEnv, ::Type{Val{:beta3}}) = env.beta3.val

getindex(env::TulipEnv, ::Type{Val{:beta4}}) = env.beta4.val

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