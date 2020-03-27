using QPSReader

# TODO: user-facing API in Julia
# Other APIs should wrap this one

# TODO: Model creation and modification
function load_problem!(m::Model{Tv}, fname::String) where{Tv}
    Base.empty!(m)

    dat = with_logger(Logging.NullLogger()) do
        readqps(fname, mpsformat=:free)
    end

    # TODO: avoid allocations when Tv is Float64
    load_problem!(m.pbdata,
        dat.name,
        true, Tv.(dat.c), Tv(dat.c0),
        sparse(dat.arows, dat.acols, Tv.(dat.avals), dat.ncon, dat.nvar),
        Tv.(dat.lcon), Tv.(dat.ucon),
        Tv.(dat.lvar), Tv.(dat.uvar),
        dat.connames, dat.varnames
    )

    return m
end


# TODO: Set/get parameters

# TODO: Query solution value and attributes

"""
    optimize!(model::Model{Tv})

Solve the optimization problem.
"""
function optimize!(model::Model{Tv}) where{Tv}
    # Instantiate solver
    model.solver = HSDSolver{Tv}(model.params, model.pbdata)

    # Solve the problem
    if isnothing(model.solver)
        error("No IPM solver is attached to the model")
        return nothing
    else
        # TODO: add a try-catch for error handling
        optimize!(model.solver, model.params)
    end

    # TODO: solution post-crush

    # Done.
    return nothing
end