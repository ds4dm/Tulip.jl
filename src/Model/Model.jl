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

    pbdata_raw::Union{Nothing, ProblemData{Tv}}   # Raw data
    pbdata_std::Union{Nothing, StandardForm{Tv}}  # Standard form

    # TODO: add the following fields:
    # * IPMSolver
    # * status
    # * ...
    Model_{Tv}() where{Tv<:Real} = new{Tv}("", ProblemData{Tv}(), nothing)

    Model_(name, pb::ProblemData{Tv}) where{Tv<:Real} = new{Tv}(name, pb, nothing)
end


include("./api.jl")