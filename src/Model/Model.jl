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
    env::Env{Tv}  # Parameters

    pbdata_raw::Union{Nothing, ProblemData{Tv}}   # Raw data
    pbdata_std::Union{Nothing, StandardForm{Tv}}  # Standard form

    # TODO: add the following fields:
    # * IPMSolver
    # * status
    # * ...
    Model{Tv}() where{Tv<:Real} = new{Tv}("", Env{Tv}(), ProblemData{Tv}(), nothing)
end


include("API/api.jl")

function optimize!(m::Model{Tv}) where{Tv<:Real}

    # convert to standard form
    sf = convert_to_standard_form(m.env.matrix_type, m.pbdata_raw)
    m.pbdata_std = sf

    # Instantiate HSD solver
    hsd = HSDSolver{Tv}(m.pbdata_std)

    # Solve problem
    optimize!(hsd, m.env)
    return hsd.solver_status
end