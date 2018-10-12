using Printf

function optimize_!(model::Model)

    # Convert to standard form
    prepross!(model)

    # Solve
    solve_ipm!(model)

end

"""
    solve_ipm!(model)

Solve the optimization problem using an interior-point algorithm.
"""
function solve_ipm!(model::Model)

    # Initialization
    # TODO: put this in a generic function
    tstart = time()
    model.time_total = 0.0
    model.num_bar_iter = 0

    # IPM log
    # TODO: put this in a generic function
    if model.env.verbose.val == 1
        println(" Itn    Primal Obj      Dual Obj        PFeas    DFeas    GFeas     Mu       Time")
    end

    # Symbolic Cholesky
    F = symbolic_cholesky(model.A)  # Symbolic factorization

    # Starting point
    # TODO: algo-specific function to generate starting point
    model.x = ones(model.num_var)
    model.w = ones(model.n_var_ub)
    model.y = zeros(model.num_constr)
    model.s = ones(model.num_var)
    model.z = ones(model.n_var_ub)
    model.t = Ref(1.0)
    model.k = Ref(1.0)
    model.μ = Ref(1.0)

    model.rp = Inf*ones(model.num_constr)
    model.ru = Inf*ones(model.n_var_ub)
    model.rd = Inf*ones(model.num_var)
    model.rg = Ref(Inf)

    # Other working memory
    θ = zeros(model.num_var)
    θwz = zeros(model.num_var_ub)

    # Main IPM loop
    # TODO (?): put this in a function (that would perform one IP iteration?)
    while (
        model.num_bar_iter < model.env.barrier_iter_max.val
        && model.time_total < model.env.time_limit.val
    )
        # ======================= #
        #     I. Book-keeping     #    
        # ======================= #

        # I.A - Compute basic info at current point
        compute_basic_info!(
            model,
            model.A,
            model.b, model.c, model.uind, model.uval,
            model.x, model.w, model.y, model.s, model.z, model.t, model.k,
            θ, θwz,
            model.rp, model.ru, model.rd, model.rg, model.μ,
            model.primbnd, model.dualbnd
        )

        # I.B - Log
        # TODO: generic display function
        model.time_total = time() - tstart
        display_ipm_log(model)

        # I.C - Check stopping criterion
        model.num_bar_iter += 1
        model.time_total = time() - tstart
        check_stopping_criterion!(model) && break

        
        # ======================== #
        #     II. Computations     #    
        # ======================== #

        # II.A - Compute Cholesky factorization
        factor_normaleq!(model.A, θ, F)
        
        # II.B - Make step
        make_step!(
            Val(model.env.algo.val),    
            model, model.env,
            model.A, model.F, model.b, model.c, model.uind, model.uval,
            θ, θwz,
            model.x, model.w, model.y, model.s, model.z, model.t, model.k
        )

    end

    return model.sln_status

end

"""
    compute_basic_info(
        A, b, c, uind, uval,
        x, w, y, s, z, t, k,
        rp, ru, rd, rg,
        μ, primbnd, dualbnd
    )

Compute residuals and basic info for HSD algorithm.

# Arguments
- `numvar`: Number of original variables
- `numvarub`: Number of upper-bounded variables
- `numcon`: Number of constraints
- `A`: Constraints' matrix
- `b`: Constraints' right-hand side
- `c`: Objective term
- `uind`: Indices of upper-bounded variables
- `uval`: Values of (finite) upper bounds on variables
- `x, w, y, s, z, t, k`: Primal and dual variables
- `rp, ru, rd, rg`: Residuals
- `μ`: Current value of the centering parameter
- `primbnd, dualbnd`: Primal and dual objective bounds
"""
function compute_basic_info!(
    numvar, numvarub, numcon,
    A::AbstractMatrix,
    b::AbstractVector, c::AbstractVector, uind, uval,
    x, w, y, s, z, t::RefValue, k::RefValue,
    θ, θwz,
    rp, ru, rd, rg::RefValue,
    μ::RefValue,
    primbnd::RefValue, dualbnd::RefValue
)

    # Diagonal scaling
    θ .= s ./ x
    θwz .= z ./ w
    @views θ[uind] .+= θwz
    θ .\= 1.0

    # Centering parameter
    μ.x = (
        (dot(x, s) + dot(w, z) + t.x * k.x)
        / (numvar + numvarub + 1)
    )

    # Primal and dual bounds
    primbnd.x = dot(c, x)
    dualbnd.x = dot(b, y) - dot(uval, z)

    # Primal residual: rp = t*b - A*x
    mul!(rp, A, x)                                      # rp = A*x
    lmul!(-oneunit(eltype(rp)), rp)                     # rp = -A*x
    axpy!(t.x, b, rp)                                   # rp = t*b - A*x

    # Upper-bound residual: ru = t*u - w - x
    mul!(ru, uval, t.x)                                 # ru = t*u
    axpy!(-oneunit(eltype(ru)), w, ru)                  # ru = t*u - w
    @views axpy!(-oneunit(eltype(ru)), x[uind], ru)     # ru = t*u - w - x

    # Dual residual: rd = t*c - A'*y - s + z
    mul!(rd, transpose(A), y)                       # rd = A'*y
    lmul!(-oneunit(eltype(rd)), rd)                 # rd = -A'*y
    axpy!(t.x, c, rd)                               # rd = t*c - A'*y
    axpy!(oneunit(eltype(rd)), s, rd)               # rd = t*c - A'*y - s
    @views rd[uind] .-= z                           # rd = t*c - A'*y - s - z

    # Gap residual: rg = c*x - b*y - u*z + k
    rg.x = primbnd.x - dualbnd.x + k.x

    return nothing

end

"""
    display_ipm_log(model)

Display IPM log for current iteration.
"""
function display_ipm_log(model::Model)

    if model.env.verbose.val == 0
        return nothing
    end

    @printf("%4d", model.num_bar_iter)
    # Primal and Dual objectives
    @printf("%+18.7e", model.primal_bound / model.t.x)
    @printf("%+16.7e", model.dual_bound / model.t.x)
    # Infeasibilities
    @printf("%10.2e", max(norm(model.rp, Inf), norm(model.ru, Inf)))  # primal infeas
    @printf("%9.2e", norm(model.rd, Inf))  # dual infeas
    @printf("%9.2e", norm(model.rg.x, Inf))  # optimality gap
    # μ
    @printf("  %7.1e", model.μ.x)
    # Time
    @printf("  %.2f", model.time_total)
    print("\n")

    return nothing

end

"""
    check_stopping_criterion!(model)

"""
function check_stopping_criterion!(model::Model)

    #================================================
        Convergence status
    ================================================# 
    # The optimization is stopped if either:
    # - Current solution is optimal
    # - A primal or dual infeasibility certificate is found
    
    eps_p = norm(model.rp, Inf) / (model.t.x + norm(model.b, Inf))   
    eps_u = norm(model.ru, Inf) / (model.t.x + norm(model.uval, Inf))
    eps_d = norm(model.rd, Inf) / (model.t.x + norm(model.c, Inf))
    eps_g = (
        abs(model.primal_bound - model.dual_bound)
        / (model.t.x + abs(model.dual_bound))
    )
    
    pfeas = (
        (eps_p <= model.env.barrier_tol_pfeas.val)
        && (eps_u <= model.env.barrier_tol_pfeas.val)
    )
    dfeas = eps_d <= model.env.barrier_tol_dfeas.val
    gfeas = eps_g <= model.env.barrier_tol_conv.val

    # TODO
    if (pfeas && dfeas && gfeas)
        # Solution is Optimal
        model.sln_status = Sln_Optimal

    # elseif
    #     # Primal-Feasible

    #     # Dual-Feasible

    #     # Primal-Dual Feasible

    #     # Primal Infeasible

    #     # Dual Infeasible

    #     # Primal-Dual Infeasible
    #     0
    else
        model.sln_status = Sln_Unknown
    end


    #================================================
        User-defined limits
    ================================================#



    return false
end

include("hsd.jl")
include("mpc.jl")