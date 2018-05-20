# import Tulip.Cholesky:
#     AbstractCholeskyFactor, SimpleDenseCholeskyFactor, SimpleSparseCholeskyFactor


"""
    solve(model, tol, verbose)

Solve model m using an infeasible predictor-corrector Interior-Point algorithm.

# Arguments
- `model::Model`: the optimization model
- `tol::Float64`: numerical tolerance
- `verbose::Int`: 0 means no output, 1 displays log at each iteration
"""
function solve!(
    model::Model;
    tol::Float64 = 10.0^-8,
    verbose::Int = 0
)

    N_ITER_MAX = 25  # maximum number of IP iterations
    niter = 0  # number of IP iterations

    # TODO: pre-optimization stuff
    F = symbolic_cholesky(model.A)  # Symbolic factorization
    θ = zeros(model.sol.x)

    # compute starting point
    compute_starting_point!(model, F)

    # TODO: check stopping criterion for possible early termination

    # IPM log
    if verbose == 1
        println(" Itn      Primal Obj        Dual Obj  Prim Inf  Dual Inf  UBnd Inf\n")
    end
        # t_fact_tot = 0.0
        # t_tot = 0.0
        # t0 = time()
    # main loop
    while niter < N_ITER_MAX
        
        # println("Current iterate:")
        # println("\t", model.sol.x)
        # println("\t", model.sol.y)
        # println("\t", model.sol.s)
        # println("\t", model.sol.w)
        # println("\t", model.sol.z)
        # println()
        # I. Form and factor Newton System
        # t1 = time()
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
            # t2 = time()
            # t_fact_tot += t2 - t1
            # t_tot = t2 - t0
            # println("t_fact = ", @sprintf("%.6f", t_fact_tot),
            #     "\tt_tot = ", @sprintf("%.6f", t_tot),
            #     "\tratio = ", @sprintf("%.6f", t_fact_tot / t_tot)
            # )

        # II. Compute and take step
        compute_next_iterate!(model, F)

        # III. Book-keeping + display log
        # compute residuals
        rb = model.A * model.sol.x - model.b
        rc = model.A' * model.sol.y + model.sol.s - model.c
        for (i,j) in enumerate(model.uind)
            rc[j] -= model.sol.z[i]
        end
    
        ru = model.sol.x[model.uind] + model.sol.w - model.uval

        obj_primal = dot(model.sol.x, model.c)
        obj_dual = dot(model.b, model.sol.y) - dot(model.uval, model.sol.z)
        # println("Current Residuals:")
        # println("\t", rb)
        # println("\t", rc)
        # println("\t", ru)
        # println()

        niter += 1
        if verbose == 1
            print(@sprintf("%4d", niter))  # iteration count
            print(@sprintf("%+18.7e", obj_primal))  # primal objective
            print(@sprintf("%+16.7e", obj_dual))  # dual objective
            print(@sprintf("%10.2e", maximum(abs.(rb))))  # primal infeas
            print(@sprintf("%9.2e", maximum(abs.(rc))))  # dual infeas
            print(@sprintf("%9.2e", maximum(abs.(ru))))  # upper bound infeas
            print("\n")
        end

        # check stopping criterion
        eps_p = (norm(rb)) / (1.0 + norm(model.b))
        eps_d = (norm(rc)) / (1.0 + norm(model.c))
        eps_u = (norm(ru)) / (1.0 + norm(model.uval))
        eps_g = abs(obj_primal - obj_dual) / (1.0 + abs(obj_primal))

        if (eps_p < tol) && (eps_u < tol) && (eps_d < tol) && (eps_g < tol)
            model.status = :Optimal
        end
        # println(
        #     @sprintf("%9.2e", eps_p),
        #     @sprintf("%9.2e", eps_u),
        #     @sprintf("%9.2e", eps_d),
        #     @sprintf("%9.2e", eps_g)
        # )

        # check status
        if model.status == :Optimal
            if verbose == 1
                println()
                println("Optimal solution found.")
            end
            return model.status
        end

    end

    # 

    return model.status
    
end


"""
    compute_starting_point!
    Compute a starting point
"""
function compute_starting_point!(model::Model, F::Factorization)
    # warn("TODO: starting point implementation")
    return model.sol
end

function compute_next_iterate!(model::Model, F::Factorization)

    (x, y, s, w, z) = (model.sol.x, model.sol.y, model.sol.s, model.sol.w, model.sol.z)
    (m, n, p) = model.nconstr, model.nvars, size(model.uind, 1)

    d_aff = copy(model.sol)
    d_cc = copy(model.sol)

    # compute residuals
    μ = (
        (dot(x, s) + dot(w, z))
        / (n + p)
    )
    rb = model.A * x - model.b
    rc = model.A' * y + s - model.c
    for (i,j) in enumerate(model.uind)
        rc[j] -= z[i]
    end

    ru = x[model.uind] + w - model.uval
    rxs = x .* s
    rwz = w .* z

    θ = x ./ s
    for (i, j) in enumerate(model.uind)
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end

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
    model.sol.x += α_p * d.x
    model.sol.y += α_d * d.y
    model.sol.s += α_d * d.s
    model.sol.w += α_p * d.w
    model.sol.z += α_d * d.z

    # println("Predictor direction:")
    # println("\td_aff.x = ", d_aff.x)
    # println("\td_aff.y = ", d_aff.y)
    # println("\td_aff.s = ", d_aff.s)
    # println("\td_aff.w = ", d_aff.w)
    # println("\td_aff.z = ", d_aff.z)
    # println()
    # println("Corrector direction:")
    # println("\t d_cc.x = ", d_cc.x)
    # println("\t d_cc.y = ", d_cc.y)
    # println("\t d_cc.s = ", d_cc.s)
    # println("\t d_cc.w = ", d_cc.w)
    # println("\t d_cc.z = ", d_cc.z)
    # println()
    # println("Final direction:")
    # println("\t    d.x = ", d.x)
    # println("\t    d.y = ", d.y)
    # println("\t    d.s = ", d.s)
    # println("\t    d.w = ", d.w)
    # println("\t    d.z = ", d.z)
    # println()

    # println("Steps:",
    #     @sprintf("%9.6f", α_pa),
    #     @sprintf("%9.6f", α_da),
    #     @sprintf("%9.6f", α_p),
    #     @sprintf("%9.6f", α_d)
    # )
    # println()
    
    return model.sol
end

"""
    symbolic_cholesky
    Compute Cholesky factorization of A*A'
"""
function symbolic_cholesky(A::AbstractMatrix{T}) where {T<:Real}

    F = cholfact(Symmetric(A*A'))
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
    for (i, j) in enumerate(uind)
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

    d.x = θ .* (A' * d.y - ξ_tmp)

    d.z = (Λ.z .* (-ξ_u + d.x[uind]) + ξ_wz) ./ Λ.w
    d.s = (ξ_xs - Λ.s .* d.x) ./ Λ.x
    d.w = (ξ_wz - Λ.w .* d.z) ./ Λ.z

    # check if system is solved correctly
    rb = ξ_b - A*d.x
    rc = ξ_c - (A'*d.y + d.s)
    rc[uind] += d.z
    ru = ξ_u - (d.x[uind] + d.w)
    rxs = ξ_xs - (Λ.s .* d.x + Λ.x .* d.s)
    rwz = ξ_wz - (Λ.z .* d.w + Λ.w .* d.z)

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

function compute_stepsize(t::PrimalDualPoint{T}, d::PrimalDualPoint{T}; damp=1.0) where T<:Real
    (ap, ad) = compute_stepsize(t.x, t.w, t.s, t.z, d.x, d.w, d.s, d.z, damp=damp)
    return (ap, ad)
end


function solve_newton_bis!(
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
    (m, n) = size(A)
    p = size(uind, 1)
    Φ = vcat(
        hcat(A, spzeros(m, p), spzeros(m, m), spzeros(m, p), spzeros(m, n)),
        hcat(speye(p), speye(p), spzeros(p, m), spzeros(p, p), spzeros(p, n)),
        hcat(spzeros(n, n), spzeros(n, p), A', -spdiagm(sparsevec(uind, ones(p), n)), speye(n)),
        hcat(spdiagm(Λ.s), spzeros(n, p), spzeros(n, m), spzeros(n, p), spdiagm(Λ.x)),
        hcat(spzeros(p, n), spdiagm(Λ.z), spzeros(p, m), spdiagm(Λ.w), spzeros(p, n))
    )

    ξ = vcat(ξ_b, ξ_u, ξ_c, ξ_xs, ξ_wz...);
    
    d_ = Φ \ ξ;
    
    d.x = d_[1:n]
    d.y = d_[(2n+1):(2n+m)]
    d.s = d_[(2n+m+p+1):end]
    d.w = d_[(n+1):(2n)]
    d.z = d_[(2n+m+1):(2n+m+p)]


    # check if system is solved correctly
    rb = ξ_b - A*d.x
    rc = ξ_c - (A'*d.y + d.s)
    rc[uind] += d.z
    ru = ξ_u - (d.x[uind] + d.w)
    rxs = ξ_xs - (Λ.s .* d.x + Λ.x .* d.s)
    rwz = ξ_wz - (Λ.z .* d.w + Λ.w .* d.z)

    println("Residuals\t(naive)")
    println("||rb||   \t", @sprintf("%.6e", maximum(abs.(rb))))
    println("||rc||   \t", @sprintf("%.6e", maximum(abs.(rc))))
    println("||ru||   \t", @sprintf("%.6e", maximum(abs.(ru))))
    println("||rxs||  \t", @sprintf("%.6e", maximum(abs.(rxs))))
    println("||rwz||  \t", @sprintf("%.6e", maximum(abs.(rwz))))
    println()

    return d
end