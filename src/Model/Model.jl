import DataStructures:
    OrderedSet, OrderedDict

"""
    ObjSense

"""
@enum ObjSense begin
    TLP_MIN  # Minimize
    TLP_MAX  # Maximize
end

import Base.*
*(s::ObjSense, x::Real) = (s == TLP_MIN) ? x : -x

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

function empty!(m::Model{Tv}) where{Tv<:Real}
    # Empty model
    m.pbdata_raw = ProblemData{Tv}()
    m.pbdata_std = nothing
    m.solver = nothing

    return m
end

function is_empty(m::Model)
    res = (
        length(m.pbdata_raw.vars) == 0
        && length(m.pbdata_raw.constrs) == 0
        && length(m.pbdata_raw.coeffs) == 0
    )
    
    return res
end

include("API/api.jl")

function optimize!(m::Model{Tv}) where{Tv<:Real}

    # Parameters
    m.env.threads > 0 || error("Invalid thread count: $(m.env.threads).")
    BLAS.set_num_threads(m.env.threads)

    # Convert to standard form
    # TODO: only re-compute what is necessary
    sf = convert_to_standard_form(m.env.matrix_type, m.pbdata_raw)
    m.pbdata_std = sf

    # Instantiate HSD solver
    # TODO: only re-compute what is necessary`
    hsd = HSDSolver{Tv}(
        m.pbdata_std.ncon, m.pbdata_std.nvar, m.pbdata_std.nupb,
        m.pbdata_std.A, m.pbdata_std.b, m.pbdata_std.c, m.pbdata_std.c0,
        m.pbdata_std.uind, m.pbdata_std.uval
    )
    m.solver = hsd

    # Solve problem
    optimize!(m.solver, m.env)

    # TODO: post-crush solution
    return m.solver.solver_status
end