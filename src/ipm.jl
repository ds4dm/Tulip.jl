import Base.LinAlg:
    A_mul_B!, At_mul_B, A_ldiv_B!

import Tulip:
    Model, PrimalDualPoint


"""
    optimize(model, tol, verbose)

Solve model m using an infeasible predictor-corrector Interior-Point algorithm.

# Arguments
- `model::Model`: the optimization model
"""
function optimize!(model::Model)

    niter = 0  # number of IP iterations

    # TODO: pre-optimization stuff
    F = symbolic_cholesky(model.A)  # Symbolic factorization
    θ = zeros(model.sol.x)

    # X = Array{PrimalDualPoint, 1}()
    # TODO: check stopping criterion for possible early termination

    # compute starting point
    # TODO: decide which starting point
    compute_starting_point!(
        model.A,
        F,
        model.sol.x,
        model.sol.w,
        model.sol.y,
        model.sol.s,
        model.sol.z,
        model.b, model.c, model.uind, model.uval
    )

    # IPM log
    if model.env[:output_level] == 1
        println(" Itn      Primal Obj        Dual Obj    Prim Inf Dual Inf UBnd Inf\n")
    end

    # main loop
    # push!(X, copy(model.sol))

    while niter < model.env[:barrier_iter_max]
        
        # I. Form and factor Newton System
        compute_newton!(
            model.A,
            model.sol.x,
            model.sol.s,
            model.sol.w,
            model.sol.z,
            model.uind,
            θ,
            F
        )

        # II. Compute and take step
        compute_next_iterate!(
            model,
            model.A,
            F,
            model.sol.x,
            model.sol.w,
            model.sol.y,
            model.sol.s,
            model.sol.z,
            model.b,
            model.c,
            model.uind,
            model.uval
        )

        # III. Book-keeping + display log
        # compute residuals
        rb = model.A * model.sol.x - model.b
        rc = At_mul_B(model.A, model.sol.y) + model.sol.s - model.c
        spxpay!(-1.0, rc, model.uind, model.sol.z)
    
        ru = model.sol.x[model.uind] + model.sol.w - model.uval

        obj_primal = dot(model.sol.x, model.c)
        obj_dual = dot(model.b, model.sol.y) - dot(model.uval, model.sol.z)

        niter += 1

        eps_p = (norm(rb)) / (1.0 + norm(model.b))
        eps_d = (norm(rc)) / (1.0 + norm(model.c))
        eps_u = (norm(ru)) / (1.0 + norm(model.uval))
        eps_g = abs(obj_primal - obj_dual) / (1.0 + abs(obj_primal))
        if model.env[:output_level] == 1
            print(@sprintf("%4d", niter))  # iteration count
            print(@sprintf("%+18.7e", obj_primal))  # primal objective
            print(@sprintf("%+16.7e", obj_dual))  # dual objective
            print(@sprintf("%10.2e", norm(rb, Inf)))  # primal infeas
            print(@sprintf("%9.2e", norm(rc, Inf)))  # dual infeas
            print(@sprintf("%9.2e", norm(ru, Inf)))  # upper bound infeas
            print(@sprintf("%9.2e", abs(obj_primal - obj_dual) / (model.n_var + model.n_var_ub)))
            print("\n")
        end

        # check stopping criterion
        if (
            (eps_p < model.env[:barrier_tol_feas])
            && (eps_u < model.env[:barrier_tol_feas])
            && (eps_d < model.env[:barrier_tol_opt])
            && (eps_g < model.env[:barrier_tol_conv])
        )
            model.status = :Optimal
        end

        # check status
        if model.status == :Optimal
            if model.env[:output_level] == 1
                println()
                println("Optimal solution found.")
            end
            return model.status
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

    d_aff = copy(model.sol)
    d_cc = copy(model.sol)
    
    # compute residuals
    μ = (
        (dot(x, s) + dot(w, z))
        / (n + p)
    )

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
        model.A,
        θ,
        F,
        model.sol,
        d_aff,
        model.uind,
        -rb,
        -rc,
        -ru,
        -rxs,
        -rwz
    )

    # compute step length
    (α_pa, α_da) = compute_stepsize(model.sol, d_aff)

    # update centrality parameter
    μ_aff = (
        (
            dot(x + α_pa * d_aff.x, s + α_da * d_aff.s)
            + dot(w + α_pa * d_aff.w, z + α_da * d_aff.z)
        ) / (n + p)
    ) 

    σ = clamp((μ_aff / μ)^3, 10.0^-12, 1.0 - 10.0^-12)  # clamped for numerical stability
    # compute corrector
    solve_newton!(
        model.A,
        θ,
        F,
        model.sol,
        d_cc,
        model.uind,
        zeros(m),
        zeros(n),
        zeros(p),
        σ*μ*ones(n) - d_aff.x .* d_aff.s,
        σ*μ*ones(p) - d_aff.w .* d_aff.z
    )

    # final step size
    d = d_aff + d_cc
    (α_p, α_d) = compute_stepsize(model.sol, d, damp=0.99995)

    # take step
    copy!(x, x + α_p * d.x)
    copy!(y, y + α_d * d.y)
    copy!(s, s + α_d * d.s)
    copy!(w, w + α_p * d.w)
    copy!(z, z + α_d * d.z)

    return model.sol
end

"""
    symbolic_cholesky
    Compute Cholesky factorization of A*A'
"""
function symbolic_cholesky(A::AbstractMatrix{T}) where {T<:Real}

    F = Cholesky.cholesky(A, ones(A.n))
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
    Cholesky.cholesky!(A, θ, F)

    return θ
end


"""
    solve_newton
    Solve Newton system with the given right-hand side.
    Overwrites the input d
"""
function solve_newton!(
    A::AbstractMatrix{Ta},
    θ::AbstractVector{T1},
    F::Factorization{Ta},
    Λ::PrimalDualPoint,
    d::PrimalDualPoint,
    uind::AbstractVector{Ti},
    ξ_b::AbstractVector{T2},
    ξ_c::AbstractVector{T3},
    ξ_u::AbstractVector{T4},
    ξ_xs::AbstractVector{T5},
    ξ_wz::AbstractVector{T6},
) where {Ta<:Real, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real, T6<:Real, Ti<:Integer}

    ξ_tmp = ξ_c - (ξ_xs ./ Λ.x)
    ξ_tmp[uind] += (ξ_wz - (Λ.z .* ξ_u)) ./ Λ.w

    d.y = F \ (ξ_b + A * (θ .* ξ_tmp))

    d.x = θ .* (At_mul_B(A, d.y) - ξ_tmp)

    d.z = (Λ.z .* (-ξ_u + d.x[uind]) + ξ_wz) ./ Λ.w
    d.s = (ξ_xs - Λ.s .* d.x) ./ Λ.x
    d.w = (ξ_wz - Λ.w .* d.z) ./ Λ.z

    # # check if system is solved correctly
    # rb = ξ_b - A*d.x
    # rc = ξ_c - (A'*d.y + d.s)
    # rc[uind] += d.z
    # ru = ξ_u - (d.x[uind] + d.w)
    # rxs = ξ_xs - (Λ.s .* d.x + Λ.x .* d.s)
    # rwz = ξ_wz - (Λ.z .* d.w + Λ.w .* d.z)

    # println("Residuals\t(normal eqs)")
    # println("||rb||   \t", @sprintf("%.6e", maximum(abs.(rb))))
    # println("||rc||   \t", @sprintf("%.6e", maximum(abs.(rc))))
    # println("||ru||   \t", @sprintf("%.6e", maximum(abs.(ru))))
    # println("||rxs||  \t", @sprintf("%.6e", maximum(abs.(rxs))))
    # println("||rwz||  \t", @sprintf("%.6e", maximum(abs.(rwz))))
    # println()
    return d
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
    
    ap, ad = -1.0, -1.0
    
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

function compute_stepsize(t::PrimalDualPoint{T}, d::PrimalDualPoint{T}; damp=1.0) where {T<:Real}
    (ap, ad) = compute_stepsize(t.x, t.w, t.s, t.z, d.x, d.w, d.s, d.z, damp=damp)
    return (ap, ad)
end

function update_theta!(θ, x, s, z, w, colind)
    # only called from within the optimization, so bounds were checked before
    for i in 1:size(colind, 1)
        j = colind[i]
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end
    return nothing
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
    return nothing
end
