import DataStructures:
    OrderedSet, OrderedDict

include("./constraint.jl")
include("./variable.jl")
include("./pbdata.jl")
include("./standardform.jl")


"""
    Model

"""
mutable struct Model
    name::String  # Model name

    # pbdata::ProblemData{Tv, Ti}

    # stdform::StandardFormData{Tv, Ti}




    # Constructor
    # function Model{Tv, Ti}(s::String) where{Tv<:Real, Ti<:Integer}
    #     return new{Tv, Ti}(s)
    # end
end