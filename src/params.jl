# Default parameters
const TLP_DEFAULT_PARAM = Dict(
    :output_level => 1,
    :barrier_iter_max => 100,
    :time_limit => Inf,
    :barrier_tol_feas => 10.0^-8,
    :barrier_tol_opt => 10.0^-8,
    :barrier_tol_conv => 10.0^-8
)

mutable struct TulipEnv

    #=======================================================
        Termination criteria
    =======================================================#

    barrier_iter_max::Int       # Maximum number of barrier iterations
    time_limit::Float64         # Time limit (in seconds)


    #=======================================================
        Tolerances
    =======================================================#

    barrier_tol_feas::Float64   # Primal feasibility tolerance
    barrier_tol_opt::Float64    # Dual feasibility tolerance
    barrier_tol_conv::Float64   # Optimality gap tolerance


    #=======================================================
        Other parameters
    =======================================================#

    output_level::Int           # 0 means no output, 1 means normal
    
    # create environment with default values
    # user can over-ride these values afterwards
    function TulipEnv()
        env = new()
        for (k, v) in TLP_DEFAULT_PARAM
            Core.setfield!(env, k, v)
        end
        return env
    end

end

function Base.copy(env::TulipEnv)
    env_ = TulipEnv()

    for (k, v) in TLP_DEFAULT_PARAM
        Core.setfield!(env_, k, Core.getfield(env, k))
    end

    return env_
end
"""
    getindex(env::TulipEnv, param)

Retrieve parameter value. Raises an error if parameter does not exist.
"""
function Base.getindex(env::TulipEnv, param::Symbol)
    if haskey(TLP_DEFAULT_PARAM, param)
        Base.getindex(env, Val{param})
    else
        error("Parameter $param does not exist")
    end
end
Base.getindex(env::TulipEnv, param::String) = Base.getindex(env, Symbol(param))


Base.getindex(env::TulipEnv, ::Type{Val{:output_level}}) = 
    copy(Core.getfield(env, :output_level))
Base.getindex(env::TulipEnv, ::Type{Val{:barrier_iter_max}}) = 
    copy(Core.getfield(env, :barrier_iter_max))
Base.getindex(env::TulipEnv, ::Type{Val{:time_limit}}) = 
    copy(Core.getfield(env, :time_limit))
Base.getindex(env::TulipEnv, ::Type{Val{:barrier_tol_feas}}) = 
    copy(Core.getfield(env, :barrier_tol_feas))
Base.getindex(env::TulipEnv, ::Type{Val{:barrier_tol_opt}}) = 
    copy(Core.getfield(env, :barrier_tol_opt))
Base.getindex(env::TulipEnv, ::Type{Val{:barrier_tol_conv}}) = 
    copy(Core.getfield(env, :barrier_tol_conv))

"""
    setindex!(env, param, v)

Set value of given parameter to v. Raises an error if parameter does not exist,
or if incorrect value.
"""
function Base.setindex!(env::TulipEnv, v, param::Symbol)

    if haskey(TLP_DEFAULT_PARAM, param)
        Base.setindex!(env, v, Val{param})
    else
        error("Parameter $param does not exist")
    end
    
    return nothing
end
Base.setindex!(env::TulipEnv, v, param::String) = setindex!(env, v, Symbol(param))




function Base.setindex!(env::TulipEnv, v, param)
    error("Function _setparam! not implemented for $param")
    return nothing
end

function Base.setindex!(env::TulipEnv, v::Real, ::Type{Val{:output_level}})

    @assert v >= 0
    env.output_level = floor(Int, v)
    return nothing
end

function Base.setindex!(env::TulipEnv, v, ::Type{Val{:barrier_iter_max}})

    @assert v >= 0
    env.barrier_iter_max = floor(Int, v)
    return nothing
end

function Base.setindex!(env::TulipEnv, v, ::Type{Val{:time_limit}})

    @assert v >= 0
    env.time_limit = Float64(v)
    return nothing
end

function Base.setindex!(env::TulipEnv, v, ::Type{Val{:barrier_tol_feas}})
    @assert 0.0 <= v <= 1.0
    env.barrier_tol_feas = Float64(v)
    return nothing
end

function Base.setindex!(env::TulipEnv, v, ::Type{Val{:barrier_tol_opt}})
    @assert 0.0 <= v <= 1.0
    env.barrier_tol_opt = Float64(v)
    return nothing
end

function Base.setindex!(env::TulipEnv, v, ::Type{Val{:barrier_tol_conv}})
    @assert 0.0 <= v <= 1.0
    env.barrier_tol_conv = Float64(v)
    return nothing
end