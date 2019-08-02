import DataStructures:
    OrderedSet, OrderedDict

include("./constraint.jl")
include("./variable.jl")
include("./pbdata.jl")


"""
    StandardFormData{Tv, Ti, Ta}

Problem data converted in standard form.
"""
mutable struct StandardFormData{Tv<:Real, Ti<:Integer}
    A::AbstractMatrix{Tv}  # Constraint matrix
    b::Vector{Tv}          # Right-hand side
    c::Vector{Tv}          # Objective
    uind::Vector{Ti}       # Indices of upper-bounded variables
    uval::Vector{Tv}       # Finite upper bounds on variables

    # TODO: add optimization sense
    # TODO: add row and column scalings
    # TODO: add starting points
end

# TODO:
# Add getters and setters
# Add some useful constructors
# Add copy and convert functions (to change numerical precision)


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