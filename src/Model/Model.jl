import OrderedCollections:
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
 
"""
    loadprolem!(m, filename)

Read problem from file `filename` and load it into `m`.
"""
function loadproblem!(m::Model{Tv}, filename::String) where{Tv<:Real}
    empty!(m)

    pb = ProblemData{Tv}()

    dat = readmps(filename)

    # Create rows
    conidx = ConstrId[]
    sizehint!(conidx, dat.ncon)
    for (i, (cname, (bt, lb, ub))) in enumerate(zip(dat.connames, dat.conbounds))
        cidx = new_constraint_index!(pb)
        constr = LinearConstraint{Tv}(cidx, cname, lb, ub)
        add_constraint!(pb, constr)
        push!(conidx, cidx)
    end

    # Create variables
    varidx = VarId[]
    sizehint!(varidx, dat.nvar)
    for (j, (vname, c, (bt, lb, ub))) in enumerate(zip(dat.varnames, dat.c, dat.varbounds))
        vidx = new_variable_index!(pb)
        var = Variable{Tv}(vidx, vname, c, lb, ub)
        add_variable!(pb, var)
        push!(varidx, vidx)
    end

    # Add coefficients
    for (i, j, v) in zip(dat.aI, dat.aJ, dat.aV)
        set_coeff!(pb, varidx[j], conidx[i], v)
    end

    # Set objective sense and offset
    if dat.objsense == :Min
        pb.obj_sense = TLP_MIN
    elseif dat.objsense == :Max
        pb.obj_sense = TLP_MAX
    else
        error("Unknown objective sense: $(dat.objsense)")
    end
    pb.obj_const = dat.c0
    
    m.name = dat.name
    m.pbdata_raw = pb

    return m
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
        m.env,
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