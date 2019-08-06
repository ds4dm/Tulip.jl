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
mutable struct Model_{Tv<:Real}
    name::String  # Model name

    # params  # Parameters

    pbdata_raw::ProblemData{Tv}   # Raw data
    # pbdata_std::StandardForm{Tv}  # Standard form

    # TODO: add the following fields:
    # * IPMSolver
    # * status
    # * ...
end

include("./api.jl")