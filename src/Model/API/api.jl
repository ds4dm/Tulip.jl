include("./vars.jl")
include("./constrs.jl")


"""
    get_value(m::Model, x::VarId)

Return value of variable `x` in current solution.
"""
function get_value(m::Model{Tv}, x::VarId) where{Tv<:Real}
    #=
        TODO
        Type inferrence breaks because we use `m.solver.pt`.
        Instead, define a function get_primal_value for IPM solvers.
    =#
    # TODO: check solver status and throw error if need be.

    # Get variable bounds in original formulation
    bt, lb, ub = get_var_bounds(m, x)

    # Index of `x` in standard form model
    j = m.pbdata_std.var2idx[x]
    x_ = m.solver.pt.x[j] / m.solver.pt.t  # Value in HSD model

    # Post-crush the solution
    # TODO: instead, perform a systematic post-crush after optimization is done
    if bt == TLP_LO || bt == TLP_FX || bt == TLP_RG
        return x_ + lb
    elseif bt == TLP_UP
        return -(x_ - ub)
    elseif bt == TLP_FR
        # Also get value of negative split
        x_m = m.solver.pt.x[j+1] / m.solver.pt.t
        return x_ - x_m
    else
        error("Variable bound $bt not supported in post-crush.")
    end
end
