include("./vars.jl")
include("./constrs.jl")


"""
    get_value(m::Model, j::VarId)

Return value of variable `x_j` in current solution.
"""
function get_value(m::Model{Tv}, j::VarId) where{Tv<:Real}
    #=
        TODO
        Type inferrence breaks because we use `m.solver.pt`.
        Instead, define a function get_primal_value for IPM solvers.
    =#
    # TODO: check solver status and throw error if need be.

    # Get variable bounds in original formulation
    bt, lb, ub = get_var_bounds(m, j)

    j_ = m.pbdata_std.var2idx[j] # Index of `x_j` in standard form model
    x = m.solver.pt.x[j_] / m.solver.pt.t  # Value in HSD model

    # Post-crush the solution
    # TODO: instead, perform a systematic post-crush after optimization is done
    if bt == TLP_LO || bt == TLP_FX || bt == TLP_RG
        return x + lb
    elseif bt == TLP_UP
        return -(x - ub)
    elseif bt == TLP_FR
        # Also get value of negative split
        x_m = m.solver.pt.x[j_+1] / m.solver.pt.t
        return x - x_m
    else
        error("Variable bound $bt not supported in post-crush.")
    end
end


"""
    get_dual_value(m::Model, i::ConstrId)

"""
function get_dual_value(m::Model{Tv}, i::ConstrId) where{Tv<:Real}
    #=
        TODO
        Type inferrence breaks because we use `m.solver.pt`.
        Instead, define a function get_dual_value for IPM solvers.
    =#
    # TODO: perform post-crush after optimization systematically
    # TODO: check solver status and throw error if need be.

    # Get index of constraint in standard form
    i_ = m.pbdata_std.con2idx[i]

    return m.solver.pt.y[i_] / m.solver.pt.t
end