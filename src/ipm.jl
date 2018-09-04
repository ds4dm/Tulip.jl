import Base.LinAlg:
    A_mul_B!, At_mul_B, A_ldiv_B!

"""
    optimize!(model)

Run the optimizer.
"""
function optimize!(model::Model)

    if model.env[:algo] == 0
        return solve_mpc!(model)
    else
        return solve_hsd!(model)
    end
end

"""
    solve_hsd(model)

Solve model using a homogeneous self-dual formulation.
"""
function solve_hsd!(model::Model)

    # Initialization
    tstart = time()
    model.runtime = 0.0
    model.numbarrieriter = 0  # number of IP iterations

    # Allocate memory for Cholesky factor
    model.F = symbolic_cholesky(model.A)  # Symbolic factorization

    # Starting point
    # TODO: enable warm-start
    model.x = ones(model.n_var)
    model.w = ones(model.n_var_ub)
    model.y = zeros(model.n_constr)
    model.s = ones(model.n_var)
    model.z = ones(model.n_var_ub)
    model.t = 1.0
    model.k = 1.0

    model.rp = Inf*ones(model.n_constr)
    model.ru = Inf*ones(model.n_var_ub)
    model.rd = Inf*ones(model.n_var)
    model.rg = Inf

    # IPM log
    if model.env[:output_level] == 1
        println(" Itn    Primal Obj      Dual Obj        PFeas    DFeas    GFeas     Mu       Time")
    end

    # main IPM loop
    while (
        model.numbarrieriter < model.env[:barrier_iter_max]
        && model.runtime < model.env[:time_limit]
    )
        # I.A - Compute residuals
        model.μ = (
            (dot(model.x, model.s) + dot(model.w, model.z) + model.t * model.k)
            / (model.n_var + model.n_var_ub + 1)
        )
        model.primal_bound = dot(model.c, model.x)
        model.dual_bound = dot(model.b, model.y) - dot(model.uval, model.z)
        compute_residuals_hsd!(
            model, model.A, model.b, model.c, model.uind, model.uval,
            model.x, model.w, model.y, model.s, model.z, model.t, model.k,
            model.rp, model.ru, model.rd, model.rg
        )
        
        # I.B - Log
        model.runtime = time() - tstart
        if model.env[:output_level] == 1
            # Iteration count
            @printf("%4d", model.numbarrieriter)
            # Primal and Dual objectives
            @printf("%+18.7e", model.primal_bound / model.t)
            @printf("%+16.7e", model.dual_bound / model.t)
            # Infeasibilities
            @printf("%10.2e", max(norm(model.rp, Inf), norm(model.ru, Inf)))  # primal infeas
            @printf("%9.2e", norm(model.rd, Inf))  # dual infeas
            @printf("%9.2e", norm(model.rg, Inf))  # optimality gap
            # μ
            @printf("  %7.1e", model.μ)
            @printf("  %.2f", model.runtime)
            print("\n")
        end

        # I.C - Stopping criteria
        if (
            (norm(model.rp, Inf) / (model.t + norm(model.b, Inf)) < model.env[:barrier_tol_feas])
            && (norm(model.ru, Inf) / (model.t + norm(model.uval, Inf)) < model.env[:barrier_tol_feas])
            && (((norm(model.rd, Inf)) / (model.t + norm(model.c, Inf))) < model.env[:barrier_tol_opt]) 
            && ((abs(model.primal_bound - model.dual_bound) / (model.t + abs(model.dual_bound))) < model.env[:barrier_tol_conv])
        )
            # optimal solution found
            if model.env[:output_level] == 1
                println("\nOptimal solution found.")
            end
            model.status = :Optimal
            return model.status
        end
        
        if (model.μ < model.env[:barrier_tol_feas]) && ((model.t / model.k) < model.env[:barrier_tol_feas])
            # infeasibility detected
            if model.env[:output_level] == 1
                println("\nInfeasibility detected.")
            end
            model.status = :Infeasible
            return model.status
        end

        # II.
        # II.A - Compute Cholesky factorization
        θ_wz = model.z ./ model.w
        θ = model.s ./ model.x
        aUtxpy!(1.0, model.uind, θ_wz, θ)
        θ .\= 1.0
        LinearAlgebra.cholesky!(model.A, θ, model.F)

        # II.B - Compute search direction
        dx, dw, dy, ds, dz, dt, dk = compute_direction_hsd(
            model, model.A, model.F, model.b, model.c, model.uval, model.uind,
            θ, θ_wz, model.μ,
            model.x, model.w, model.y, model.s, model.z, model.t, model.k,
            model.rp, model.ru, model.rd, model.rg
        )

        # II.C - Make step
        make_step_hsd!(
            model,
            model.x, model.w, model.y, model.s, model.z, model.t, model.k,
            dx, dw, dy, ds, dz, dt, dk
        )

        model.numbarrieriter += 1

    end

    # END
    return model.status
    
end

"""
    solve_mpc(model)

Solve model using Mehrothra's predictor-corrector algorithm.
"""
function solve_mpc!(model::Model)

    # Initialization
    tstart = time()
    model.runtime = 0.0
    model.numbarrieriter = 0  # number of IP iterations

    # Initialize iterates
    model.F = symbolic_cholesky(model.A)  # Symbolic factorization
    θ = zeros(model.x)
    model.x = Vector{Float64}(model.n_var)
    model.w = Vector{Float64}(model.n_var_ub)
    model.y = Vector{Float64}(model.n_constr)
    model.s = Vector{Float64}(model.n_var)
    model.z = Vector{Float64}(model.n_var_ub)

    # compute starting point
    # TODO: decide which starting point
    compute_starting_point!(
        model.A,
        model.F,
        model.x,
        model.w,
        model.y,
        model.s,
        model.z,
        model.b,
        model.c,
        model.uind,
        model.uval
    )

    # IPM log
    if model.env[:output_level] == 1
        println(" Itn    Primal Obj      Dual Obj        Prim Inf Dual Inf UBnd Inf")
    end

    # main IPM loop
    while (
        model.numbarrieriter < model.env[:barrier_iter_max]
        && model.status != :Optimal
        && model.runtime < model.env[:time_limit]
    )

        # I. Form and factor Newton System
        compute_newton!(
            model.A,
            model.x,
            model.s,
            model.w,
            model.z,
            model.uind,
            θ,
            model.F
        )

        # II. Compute and take step
        compute_next_iterate!(
            model,
            model.A,
            model.F,
            model.x,
            model.w,
            model.y,
            model.s,
            model.z,
            model.b,
            model.c,
            model.uind,
            model.uval
        )

        # III. Book-keeping + display log
        # compute residuals
        rp = model.A * model.x - model.b
        rd = At_mul_B(model.A, model.y) + model.s - model.c
        aUtxpy!(-1.0, model.uind, model.z, rd)

        ru = model.x[model.uind] + model.w - model.uval

        obj_primal = dot(model.x, model.c)
        obj_dual = dot(model.b, model.y) - dot(model.uval, model.z)

        model.numbarrieriter += 1

        eps_p = (norm(rp)) / (1.0 + norm(model.b))
        eps_d = (norm(rd)) / (1.0 + norm(model.c))
        eps_u = (norm(ru)) / (1.0 + norm(model.uval))
        eps_g = abs(obj_primal - obj_dual) / (1.0 + abs(obj_primal))


        # check stopping criterion
        if (
            (eps_p < model.env[:barrier_tol_feas])
            && (eps_u < model.env[:barrier_tol_feas])
            && (eps_d < model.env[:barrier_tol_opt])
            && (eps_g < model.env[:barrier_tol_conv])
        )
            model.status = :Optimal
        end

        # Log
        model.runtime = time() - tstart

        if model.env[:output_level] == 1
            # Iteration count
            print(@sprintf("%4d", model.numbarrieriter))
            # Primal and Dual objectives
            print(@sprintf("%+18.7e", obj_primal))
            print(@sprintf("%+16.7e", obj_dual))
            # Infeasibilities
            print(@sprintf("%10.2e", norm(rp, Inf)))  # primal infeas
            print(@sprintf("%9.2e", norm(rd, Inf)))  # dual infeas
            print(@sprintf("%9.2e", norm(ru, Inf)))  # upper bound infeas
            # μ
            print(@sprintf("  %8.2e", abs(obj_primal - obj_dual) / (model.n_var + model.n_var_ub)))
            print(@sprintf("  %.2f", model.runtime))
            print("\n")
        end

    end

    # END
    return model.status
    
end


"""
    compute_residuals_hsd!

Compute residuals at the current point
"""
function compute_residuals_hsd!(
    model, A, b, c, uind, uval,
    x, w, y, s, z, t, k,
    rp, ru, rd, rg
)
    # primal residual
    rp .= b*t - A*x

    # upper-bound residual
    ru .= uval * t - w
    aUxpy!(-1.0, x, uind, ru)

    # dual residual
    rd .= c*t - At_mul_B(A, y) - s
    aUtxpy!(1.0, uind, z, rd)

    # gap residual
    model.rg = model.primal_bound - model.dual_bound + k

    return nothing
    
end

function compute_direction_hsd(
    model, A, F, b, c, uval, uind,
    θ, θ_wz, μ,
    x, w, y, s, z, t, k,
    rp, ru, rd, rg
)
    # compute predictor direction
    dx_a, dw_a, dy_a, ds_a, dz_a, dt_a, dk_a = solve_newton_hsd(
        A, F, b, c, uval, uind,
        θ, θ_wz,
        x, w, y, s, z, t, k,
        rp, ru, rd, rg, -x .* s, -w .* z, -t*k
    )

    
    a = compute_max_step_size(
        x, w, y, s, z, t, k,
        dx_a, dw_a, dy_a, ds_a, dz_a, dt_a, dk_a
    )
    
    μ_a = (
        dot(x+a*dx_a, s+a*ds_a)
        + dot(w + a * dw_a, z + a*dz_a)
        + (t+a*dt_a)*(k+a*dk_a)
    ) / (model.n_var + model.n_var_ub + 1)

    γ = (1-a)^2 * min(1-a, 0.1)  # TODO: replace 0.1 by parameter β1
    η = 1.0 - γ

    # compute corrector
    dx, dw, dy, ds, dz, dt, dk = solve_newton_hsd(
        A, F, b, c, uval, uind,
        θ, θ_wz,
        x, w, y, s, z, t, k,
        η*rp, η*ru, η*rd, η*rg,
        -x .* s - dx_a .* ds_a + γ * μ * ones(model.n_var),
        -w .* z - dw_a .* dz_a + γ * μ * ones(model.n_var_ub),
        -t*k - dt_a * dk_a + γ*μ
    )

    return dx, dw, dy, ds, dz, dt, dk
end

"""
    compute_starting_point!

Compute a starting point

# Arguments
-`model::Model`
-`F::Factorization`: Cholesky factor of A*A', where A is the constraint matrix
    of the optimization problem.
"""
function compute_starting_point!(
    A::AbstractMatrix{Tv},
    F::Factorization{Tv},
    x::AbstractVector{Tv},
    w::AbstractVector{Tv},
    y::AbstractVector{Tv},
    s::AbstractVector{Tv},
    z::AbstractVector{Tv},
    b::StridedVector{Tv},
    c::StridedVector{Tv},
    uind::StridedVector{Ti},
    uval::StridedVector{Tv}
) where{Tv<:Real, Ti<:Integer}

    (m, n) = size(A)
    p = size(uind, 1)

    rhs = - 2 * b
    u_ = zeros(n)
    aUtxpy!(1.0, uind, uval, u_)
    rhs += A* u_


    #=======================================================
        I. Compute initial points
    =======================================================#

    # Compute x0
    v = F \ rhs

    copy!(x, -0.5 * At_mul_B(A, v))
    aUtxpy!(0.5, uind, uval, x)

    # Compute w0
    @inbounds for i in 1:p
        j = uind[i]
        w[i] = uval[i] - x[j]
    end

    # Compute y0
    copy!(y, F \ (A*c))

    # Compute s0
    copy!(s, 0.5 * (At_mul_B(A, y) - c))

    # Compute z0
    @inbounds for i in 1:p
        j = uind[i]
        z[i] = - s[j]
    end


    #=======================================================
        II. Correction
    =======================================================# 

    dp = zero(Tv)
    dd = zero(Tv)

    @inbounds for i in 1:n
        tmp = -1.5 * x[i]
        if tmp > dp
            dp = tmp
        end
    end

    @inbounds for i in 1:p
        tmp = -1.5 * w[i]
        if tmp > dp
            dp = tmp
        end
    end


    @inbounds for i in 1:n
        tmp = -1.5 * s[i]
        if tmp > dd
            dd = tmp
        end
    end

    @inbounds for i in 1:p
        tmp = -1.5 * z[i]
        if tmp > dd
            dd = tmp
        end
    end

    tmp = dot(x + dp, s + dd) + dot(w + dp, z + dd)

    dp += 0.5 * tmp / (sum(s + dd) + sum(z + dd))
    dd += 0.5 * tmp / (sum(x + dp) + sum(w + dp))

    #=======================================================
        III. Apply correction
    =======================================================#

    @inbounds for i in 1:n
        x[i] += dp    
    end

    @inbounds for i in 1:n
        s[i] += dd    
    end

    @inbounds for i in 1:p
        w[i] += dp    
    end

    @inbounds for i in 1:p
        z[i] += dd    
    end

    # Done
    return nothing
end

function compute_next_iterate!(
    model,
    A::AbstractMatrix{Tv},
    F::Factorization{Tv},
    x::AbstractVector{Tv},
    w::AbstractVector{Tv},
    y::AbstractVector{Tv},
    s::AbstractVector{Tv},
    z::AbstractVector{Tv},
    b::StridedVector{Tv},
    c::StridedVector{Tv},
    uind::StridedVector{Ti},
    uval::StridedVector{Tv}
) where{Tv<:Real, Ti<:Integer}
    (m, n, p) = model.n_constr, model.n_var, model.n_var_ub

    # Affine-scaling direction
    daff_x = copy(x)
    daff_w = copy(w)
    daff_y = copy(y)
    daff_s = copy(s)
    daff_z = copy(z)

    # Corrector-Centering direction
    dcc_x = copy(x)
    dcc_w = copy(w)
    dcc_y = copy(y)
    dcc_s = copy(s)
    dcc_z = copy(z)
    
    # compute residuals
    μ = (dot(x, s) + dot(w, z)) / (n + p)

    rp = (A * x) - b
    rd = At_mul_B(A, y) + s - c
    aUtxpy!(-1.0, uind, z, rd)

    ru = x[uind] + w - uval
    rxs = x .* s
    rwz = w .* z

    θ = x ./ s
    update_theta!(θ, x, s, z, w, uind)

    # compute predictor
    solve_newton!(
        A,
        θ,
        F,
        x, w, y, s, z,
        daff_x, daff_w, daff_y, daff_s, daff_z,
        uind,
        -rp,
        -rd,
        -ru,
        -rxs,
        -rwz
    )

    # compute step length
    (α_pa, α_da) = compute_stepsize(
        x, w, s, z,
        daff_x, daff_w, daff_s, daff_z 
    )

    # update centrality parameter
    μ_aff = (dot(x + α_pa * daff_x, s + α_da * daff_s)
            +dot(w + α_pa * daff_w, z + α_da * daff_z)) / (n + p)

    σ = clamp((μ_aff / μ)^3, 10.0^-12, 1.0 - 10.0^-12)  # clamped for numerical stability
    # compute corrector
    solve_newton!(
        A,
        θ,
        F,
        x, w, y, s, z,
        dcc_x, dcc_w, dcc_y, dcc_s, dcc_z,
        uind,
        zeros(m),
        zeros(n),
        zeros(p),
        σ*μ*ones(n) - daff_x .* daff_s,
        σ*μ*ones(p) - daff_w .* daff_z
    )

    # final step size
    dx = daff_x + dcc_x
    dw = daff_w + dcc_w
    dy = daff_y + dcc_y
    ds = daff_s + dcc_s
    dz = daff_z + dcc_z

    (α_p, α_d) = compute_stepsize(
        x, w, s, z,
        dx, dw, ds, dz,
        damp=0.99995
    )

    # take step
    x .+= α_p * dx
    w .+= α_p * dw
    y .+= α_d * dy
    s .+= α_d * ds
    z .+= α_d * dz

    return nothing
end

"""
    symbolic_cholesky
    Compute Cholesky factorization of A*A'
"""
function symbolic_cholesky(A::AbstractMatrix{T}) where {T<:Real}

    F = LinearAlgebra.cholesky(A, ones(A.n))
    return F

end


"""
    compute_newton!
    Form and factorize the Newton system, using the normal equations.
"""
function compute_newton!(
        A::AbstractMatrix{Ta},
        x::AbstractVector{Tx},
        s::AbstractVector{Ts},
        w::AbstractVector{Tw},
        z::AbstractVector{Tz},
        uind::AbstractVector{Ti},
        θ::AbstractVector{T},
        F::Factorization{Ta}
    ) where {Ta<:Real, Tx<:Real, Ts<:Real, Tw<:Real, Tz<:Real, Ti<:Integer, T<:Real}

    # Compute Θ = (X^{-1} S + W^{-1} Z)^{-1}
    θ = x ./ s
    for i in 1:size(uind, 1)
        j = uind[i]
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end

    # Form the normal equations matrix and compute its factorization
    LinearAlgebra.cholesky!(A, θ, F)

    return θ
end


"""
    solve_newton
    Solve Newton system with the given right-hand side.
    Overwrites the input d
"""
function solve_newton!(
    A::Ta,
    θ::AbstractVector{Tv},
    F::Factorization{Tv},
    x::AbstractVector{Tv},
    w::AbstractVector{Tv},
    y::AbstractVector{Tv},
    s::AbstractVector{Tv},
    z::AbstractVector{Tv},
    dx::AbstractVector{Tv},
    dw::AbstractVector{Tv},
    dy::AbstractVector{Tv},
    ds::AbstractVector{Tv},
    dz::AbstractVector{Tv},
    uind::AbstractVector{Int},
    ξ_b::AbstractVector{Tv},
    ξ_c::AbstractVector{Tv},
    ξ_u::AbstractVector{Tv},
    ξ_xs::AbstractVector{Tv},
    ξ_wz::AbstractVector{Tv},
) where {Tv<:Real, Ta<:AbstractMatrix{Tv}}

    ξ_tmp = ξ_c - (ξ_xs ./ x)
    ξ_tmp[uind] += (ξ_wz - (z .* ξ_u)) ./ w

    dy .= F \ (ξ_b + A * (θ .* ξ_tmp))

    dx .= θ .* (At_mul_B(A, dy) - ξ_tmp)

    dz .= (z .* (-ξ_u + dx[uind]) + ξ_wz) ./ w
    ds .= (ξ_xs - s .* dx) ./ x
    dw .= (ξ_wz - w .* dz) ./ z

    return nothing
end

function solve_newton_hsd(
    A, F, b, c, uval, uind,
    θ, θ_wz,
    x, w, y, s, z, t, k,
    rp, ru, rd, rg, rxs, rwz, rtk
)
    # Compute γ1, γ2, γ3
    γ_x, γ_y, γ_z = solve_augmented_system_hsd(A, F, θ, θ_wz, uind, b, c, uval)
    
    # Compute δ1, δ2, δ3
    δ_x, δ_y, δ_z = solve_augmented_system_hsd(A, F, θ, θ_wz, uind, rp, rd - rxs ./ x, ru - rwz ./ z)
    
    # Compute Δτ
    dt = (
        (rg + (rtk / t) + dot(c, δ_x) - dot(b, δ_y) + dot(uval, δ_z) ) 
        / ((k / t) - dot(c, γ_x) + dot(b, γ_y) - dot(uval, γ_z))
    )
    
    # Compute Δx, Δy
    dx = δ_x + dt * γ_x
    dy = δ_y + dt * γ_y
    dz = δ_z + dt * γ_z
    
    # Compute Δs, Δκ, Δw
    dw = (rwz - w .* dz) ./ z
    ds = (rxs - s .* dx) ./ x
    dk = (rtk - k * dt) / t

    return dx, dw, dy, ds, dz, dt, dk
end

function solve_augmented_system_hsd(A, F, θ, θ_wz, uind, rp, rd, ru)
    # right-hand side
    ru_ = ru .* θ_wz  # (W^{-1}*Z)*ru
    
    rp_ = copy(rd)
    aUtxpy!(-1.0, uind, ru_, rp_)
    rp_ .*= θ
    rhs_y = rp + A * rp_
    
    # compute y
    y = F \ rhs_y
    x = At_mul_B(A, y) - rd
    aUtxpy!(1.0, uind, ru_, x)
    x .*= θ
    
    z = copy(-ru)
    aUxpy!(1.0, x, uind, z)
    z .*= θ_wz
    
    return x, y, z
end

"""
    elementary_step_size(x, dx)

Compute maximum step size `0 <= a` such that `x + a*dx >= 0`
"""
function elementary_step_size(x, dx)
    n = size(x, 1)
    n == size(dx, 1) || throw(DimensionMismatch())
    a = Inf
    @inbounds for i in 1:n
        if dx[i] < 0
            if (-x[i] / dx[i]) < a
                a = (-x[i] / dx[i])
            end
        end
    end
    return a
end

function compute_max_step_size(x, w, y, s, z, t, k, dx, dw, dy, ds, dz, dt, dk)
    
    ax = elementary_step_size(x, dx)
    aw = elementary_step_size(w, dw)
    as = elementary_step_size(s, ds)
    az = elementary_step_size(z, dz)
    
    at = dt < 0.0 ? (-t / dt) : 1.0
    ak = dk < 0.0 ? (-k / dk) : 1.0
    
    a = min(1.0, ax, aw, as, az, at, ak)
    
    return a

end

function compute_step_size_hsd(x, w, y, s, z, t, k, dx, dw, dy, ds, dz, dt, dk)
    
    ax = elementary_step_size(x, dx)
    aw = elementary_step_size(w, dw)
    as = elementary_step_size(s, ds)
    az = elementary_step_size(z, dz)
    
    at = dt < 0.0 ? (-t / dt) : Inf
    ak = dk < 0.0 ? (-k / dk) : Inf
    
    ap = 0.99995 * min(1.0, ax, aw, at, ak)
    ad = 0.99995 * min(1.0, as, az, at, ak)
    
    return ap, ad

end

function compute_stepsize(
    x::AbstractVector{T}, w::AbstractVector{T}, s::AbstractVector{T}, z::AbstractVector{T},
    dx::AbstractVector{T}, dw::AbstractVector{T}, ds::AbstractVector{T}, dz::AbstractVector{T};
    damp=1.0
) where T<:Real
    
    ax = elementary_step_size(x, dx)
    aw = elementary_step_size(w, dw)
    as = elementary_step_size(s, ds)
    az = elementary_step_size(z, dz)
    
    ap = damp * min(ax, aw, 1.0)
    ad = damp * min(as, az, 1.0)
    
    return (ap, ad)

end

function make_step_hsd!(
    model, x, w, y, s, z, t, k,
    dx, dw, dy, ds, dz, dt, dk
)

    ap, ad = compute_step_size_hsd(x, w, y, s, z, t, k, dx, dw, dy, ds, dz, dt, dk)
    a = min(ap, ad)
    x .+= a * dx
    w .+= a * dw
    y .+= a * dy
    s .+= a * ds
    z .+= a * dz
    model.t += a * dt
    model.k += a * dk

    return nothing
end

function update_theta!(θ, x, s, z, w, colind)
    # only called from within the optimization, so bounds were checked before
    for i in 1:size(colind, 1)
        j = colind[i]
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end
end

"""
    aUtxpy!(a, xval, xind, y)

Lift `x` to a `n`-dimensional space, and add the scaled result to `y`.
Equivalent to writing `y[xind] .+= a * xval`.
"""
function aUtxpy!(a, xind, xval, y)
    
    n = size(y, 1)
    p = size(xind, 1)
    p == size(xval, 1) || throw(DimensionMismatch("xind has size $p, but xval has size $(size(xval))"))
    if p == 0 || a == zero(a)
        # early return
        return y
    end
    n >= xind[end] || throw(DimensionMismatch("Too many indices: $(xind[end]) > $n"))
    
    
    
    if a == oneunit(a)
        @inbounds for (j, i) in enumerate(xind)
            y[i] += xval[j]  # i = xind[j]
        end
    else
        @inbounds for (j, i) in enumerate(xind)
            y[i] += a * xval[j]  # i = xind[j]
        end
    end
    return y
end

"""
    aUxpy!(a, x, yind, yval)

Project `x` down to a lower-dimensional space, and add the scaled result to `y`.
Equivalent to writing `yval .+= a * x[yind]`
"""
function aUxpy!(a, x, yind, yval)
    
    n = size(x, 1)
    p = size(yind, 1)
    p == size(yval, 1) || throw(DimensionMismatch())
    if p==0 || a == zero(a)
        # early return
        return yval
    end

    n >= yind[end] || throw(DimensionMismatch())

    
    
    if a == oneunit(a)
        @inbounds for (j, i) in enumerate(yind)
            yval[j] += x[i]  # i = xind[j]
        end
    else
        @inbounds for (j, i) in enumerate(yind)
            yval[j] += a * x[i]  # i = xind[j]
        end
    end
    
    return yval
end