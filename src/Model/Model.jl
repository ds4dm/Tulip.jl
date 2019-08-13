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

    solver::Union{Nothing, AbstractIPMSolver{Tv}}
    # TODO: add the following fields:
    # * IPMSolver
    # * status
    # * ...
    function Model{Tv}() where{Tv<:Real}
        m = new{Tv}()
        m.name = ""
        m.env = Env{Tv}()
        
        m.pbdata_raw = ProblemData{Tv}()
        m.pbdata_std = nothing
        m.solver = nothing
        return m
    end
end


include("API/api.jl")

function optimize!(m::Model{Tv}) where{Tv<:Real}

    if isa(m.pbdata_std, Nothing)
        # convert to standard form
        sf = convert_to_standard_form(m.env.matrix_type, m.pbdata_raw)
        m.pbdata_std = sf
    end

    if isa(m.solver, Nothing)
        # Instantiate HSD solver
        hsd = HSDSolver{Tv}(
            m.pbdata_std.ncon, m.pbdata_std.nvar, m.pbdata_std.nupb,
            m.pbdata_std.A, m.pbdata_std.b, m.pbdata_std.c,
            m.pbdata_std.uind, m.pbdata_std.uval
        )
        m.solver = hsd
    end

    # Solve problem
    # TODO: un-crunch solution after optimization
    optimize!(m.solver, m.env)
    return m.solver.solver_status
end