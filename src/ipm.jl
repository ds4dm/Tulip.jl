import Base.LinAlg:
    A_mul_B!, At_mul_B, A_ldiv_B!

"""
    optimize(model, tol, verbose)

Solve model m using an infeasible predictor-corrector Interior-Point algorithm.

# Arguments
- `model::Model`: the optimization model
"""
function optimize!(model::Model)

    tstart = time()
    model.runtime = 0.0
    model.numbarrieriter = 0  # number of IP iterations

    # TODO: pre-optimization stuff
    model.F = symbolic_cholesky(model.A)  # Symbolic factorization
    θ = zeros(model.x)

    # X = Array{PrimalDualPoint, 1}()
    # TODO: check stopping criterion for possible early termination

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
        rb = model.A * model.x - model.b
        rc = At_mul_B(model.A, model.y) + model.s - model.c
        spxpay!(-1.0, rc, model.uind, model.z)

        ru = model.x[model.uind] + model.w - model.uval

        obj_primal = dot(model.x, model.c)
        obj_dual = dot(model.b, model.y) - dot(model.uval, model.z)

        model.numbarrieriter += 1

        eps_p = (norm(rb)) / (1.0 + norm(model.b))
        eps_d = (norm(rc)) / (1.0 + norm(model.c))
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
            print(@sprintf("%10.2e", norm(rb, Inf)))  # primal infeas
            print(@sprintf("%9.2e", norm(rc, Inf)))  # dual infeas
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
    spxpay!(1.0, u_, uind, uval)
    rhs += A* u_


    #=======================================================
        I. Compute initial points
    =======================================================#

    # Compute x0
    v = F \ rhs

    copy!(x, -0.5 * At_mul_B(A, v))
    spxpay!(0.5, x, uind, uval)

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
    (m, n, p) = model.n_con, model.n_var, model.n_var_ub

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

    rb = (A * x) - b
    rc = At_mul_B(A, y) + s - c
    spxpay!(-1.0, rc, uind, z)

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
        -rb,
        -rc,
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


function compute_stepsize(
    tx::AbstractVector{T}, tw::AbstractVector{T}, ts::AbstractVector{T}, tz::AbstractVector{T},
    dx::AbstractVector{T}, dw::AbstractVector{T}, ds::AbstractVector{T}, dz::AbstractVector{T};
    damp=1.0
) where T<:Real

    n = size(tx, 1)
    p = size(tw, 1)
    n == size(ts, 1) || throw(DimensionMismatch("t.s is wrong size"))
    p == size(tz, 1) || throw(DimensionMismatch("t.z is wrong size"))
    
    n == size(dx, 1) || throw(DimensionMismatch("d.x is wrong size"))
    n == size(ds, 1) || throw(DimensionMismatch("d.s is wrong size"))
    p == size(dw, 1) || throw(DimensionMismatch("d.w is wrong size"))
    p == size(dz, 1) || throw(DimensionMismatch("d.z is wrong size"))
    
    ap = ad = -1.0
    
    @inbounds for i in 1:n
        if dx[i] < 0.0
            if (tx[i] / dx[i]) > ap
                ap = (tx[i] / dx[i])
            end
        end
    end
    
    @inbounds for i in 1:n
        if ds[i] < 0.0
            if (ts[i] / ds[i]) > ad
                ad = (ts[i] / ds[i])
            end
        end
    end
    
    @inbounds for j in 1:p
        if dw[j] < 0.0
            if (tw[j] / dw[j]) > ap
                ap = (tw[j] / dw[j])
            end
        end
    end
    
    @inbounds for j in 1:p
        if dz[j] < 0.0
            if (tz[j] / dz[j]) > ad
                ad = (tz[j] / dz[j])
            end
        end
    end
    
    ap = - damp * ap
    ad = - damp * ad

    # println("\t(ap, ad) = ", (ap, ad))
    
    return (ap, ad)

end

function update_theta!(θ, x, s, z, w, colind)
    # only called from within the optimization, so bounds were checked before
    for i in 1:size(colind, 1)
        j = colind[i]
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end
end

"""
    spxpay!(α, x, y_ind, y_val)

In-place computation of x += α * y, where y = sparse(y_ind, y_val)

# Arguments
"""
function spxpay!(α::Tv, x::AbstractVector{Tv}, y_ind::AbstractVector{Ti}, y_val::AbstractVector{Tv}) where{Ti<:Integer, Tv<:Real}
    for i in 1:size(y_ind, 1)
        j = y_ind[i]
        x[j] = x[j] + α * y_val[i]
    end
end
