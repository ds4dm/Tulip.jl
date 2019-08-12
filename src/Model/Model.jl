import DataStructures:
    OrderedSet, OrderedDict

include("./constraint.jl")
include("./variable.jl")

# Q: add AbstractFormulation{Tv} type and methods?
include("./pbdata.jl")
include("./standardform.jl")


"""
    Model{Tv}

"""
mutable struct Model{Tv<:Real}
    name::String  # Model name
    env::TulipEnv  # Parameters

    pbdata_raw::Union{Nothing, ProblemData{Tv}}   # Raw data
    pbdata_std::Union{Nothing, StandardForm{Tv}}  # Standard form

    # TODO: add the following fields:
    # * IPMSolver
    # * status
    # * ...
    Model{Tv}() where{Tv<:Real} = new{Tv}("", TulipEnv(), ProblemData{Tv}(), nothing)
end


include("API/api.jl")