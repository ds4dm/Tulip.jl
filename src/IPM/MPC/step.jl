"""
    compute_step!(ipm, params)

Compute next IP iterate for the MPC formulation.

# Arguments
- `ipm`: The MPC optimizer model
- `params`: Optimization parameters
"""
function compute_step!(mpc::MPC{T, Tv}, params::IPMOptions{T}) where{T, Tv<:AbstractVector{T}}

    # Names
    dat = mpc.dat
    pt = mpc.pt
    res = mpc.res

    m, n, p = pt.m, pt.n, pt.p

    A = dat.A
    b = dat.b
    c = dat.c

    # Compute scaling
    θl = (pt.zl ./ pt.xl) .* dat.lflag
    θu = (pt.zu ./ pt.xu) .* dat.uflag
    θinv = θl .+ θu

    # Update regularizations
    mpc.regP .= clamp(sqrt(pt.μ), sqrt(eps(T)), sqrt(sqrt(eps(T))))
    mpc.regD .= clamp(sqrt(pt.μ), sqrt(eps(T)), sqrt(sqrt(eps(T))))

    # Update factorization
    nbump = 0
    while nbump <= 3
        try
            @timeit mpc.timer "Factorization" KKT.update!(mpc.kkt, θinv, mpc.regP, mpc.regD)
            break
        catch err
            isa(err, PosDefException) || isa(err, ZeroPivotException) || rethrow(err)

            # Increase regularization
            mpc.regD .*= 100
            mpc.regP .*= 100
            nbump += 1
            @warn "Increase regularizations to $(mpc.regP[1])"
        end
    end
    # TODO: throw a custom error for numerical issues
    nbump < 3 || throw(PosDefException(0))  # factorization could not be saved

    # Search directions
    # Predictor
    Δ  = Point{T, Tv}(m, n, p, hflag=false)
    Δc = Point{T, Tv}(m, n, p, hflag=false)

    # Affine-scaling direction
    @timeit mpc.timer "Newton" solve_newton_system!(Δ, mpc,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd,
        -(pt.xl .* pt.zl) .* dat.lflag,
        -(pt.xu .* pt.zu) .* dat.uflag,
    )

    # Step length for affine-scaling direction
    α = max_step_length(pt, Δ)
    μₐ = (
        dot((pt.xl + α .* Δ.xl) .* dat.lflag, pt.zl + α .* Δ.zl)
        + dot((pt.xu + α .* Δ.xu) .* dat.uflag, pt.zu + α .* Δ.zu)
    ) / pt.p
    σ = clamp((μₐ / pt.μ)^3, sqrt(eps(T)), one(T) - sqrt(eps(T)))
    
    # Mehrotra corrector
    @timeit mpc.timer "Newton" solve_newton_system!(Δ, mpc,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd,
        (-(pt.xl .* pt.zl) .+ σ * pt.μ .- Δ.xl .* Δ.zl) .* dat.lflag,
        (-(pt.xu .* pt.zu) .+ σ * pt.μ .- Δ.xu .* Δ.zu) .* dat.uflag,
    )
    α = max_step_length(pt, Δ)

    # Extra corrections
    #=ncor = 0
    while ncor < params.BarrierCorrectionLimit && α < T(999 // 1000)
        α_ = α
        ncor += 1

        # Compute extra-corrector
        αc = compute_higher_corrector!(Δc, 
            mpc, γ, 
            hx, hy, h0,
            Δ, α_, params.BarrierCentralityOutlierThreshold
        )
        if αc > α_
            # Use corrector
            Δ.x  .= Δc.x
            Δ.xl .= Δc.xl
            Δ.xu .= Δc.xu
            Δ.y  .= Δc.y
            Δ.zl .= Δc.zl
            Δ.zu .= Δc.zu
            Δ.τ   = Δc.τ
            Δ.κ   = Δc.κ
            α = αc
        end
        
        if αc < T(11 // 10) * α_
            break
        end

        # if (1.0 - η * αc) >= 0.9*(1.0 - η * α_)
        #     # not enough improvement, step correcting
        #     break
        # end
        
    end =#

    # Update current iterate
    α *= params.BarrierStepDampFactor
    pt.x  .+= α .* Δ.x
    pt.xl .+= α .* Δ.xl
    pt.xu .+= α .* Δ.xu
    pt.y  .+= α .* Δ.y
    pt.zl .+= α .* Δ.zl
    pt.zu .+= α .* Δ.zu
    update_mu!(pt)
    
    return nothing
end


"""
    solve_newton_system!(Δ, mpc, ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk)

Solve the Newton system
```math
\\begin{bmatrix}
    A & & & R_{d} & & \\\\
    I & -I & & & & \\\\
    I & & I & & & \\\\
    -R_{p} & & & A^{T} & I & -I \\\\
    & Z_{l} & & & X_{l}\\\\
    & & Z_{u} & & & X_{u}\\\\
\\end{bmatrix}
\\begin{bmatrix}
    Δ x\\\\
    Δ x_{l}\\\\
    Δ x_{u}\\\\
    Δ y\\\\
    Δ z_{l} \\\\
    Δ z_{u}
\\end{bmatrix}
=
\\begin{bmatrix}
    ξ_p\\\\
    ξ_l\\\\
    ξ_u\\\\
    ξ_d\\\\
    ξ_{xz}^{l}\\\\
    ξ_{xz}^{u}
\\end{bmatrix}
```

# Arguments
- `Δ`: Search direction, modified
- `mpc`: The MPC optimizer
- `hx, hy, hz, h0`: Terms obtained in the preliminary augmented system solve
- `ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk`: Right-hand side vectors
"""
function solve_newton_system!(Δ::Point{T, Tv},
    mpc::MPC{T, Tv},
    # Right-hand side
    ξp::Tv, ξl::Tv, ξu::Tv, ξd::Tv, ξxzl::Tv, ξxzu::Tv
) where{T, Tv<:AbstractVector{T}}

    pt = mpc.pt
    dat = mpc.dat

    # I. Solve augmented system
    @timeit mpc.timer "ξd_"  begin
        ξd_ = copy(ξd)
        @. ξd_ += -((ξxzl + pt.zl .* ξl) ./ pt.xl) .* dat.lflag + ((ξxzu - pt.zu .* ξu) ./ pt.xu) .* dat.uflag
    end
    @timeit mpc.timer "KKT" KKT.solve!(Δ.x, Δ.y, mpc.kkt, ξp, ξd_)

    # II. Recover Δxl, Δxu
    @timeit mpc.timer "Δxl" begin
        @. Δ.xl = -ξl + Δ.x
        Δ.xl .*= dat.lflag
    end
    @timeit mpc.timer "Δxu" begin
        @. Δ.xu =  ξu - Δ.x
        Δ.xu .*= dat.uflag
    end

    # III. Recover Δzl, Δzu
    @timeit mpc.timer "Δzl" @. Δ.zl = ((ξxzl - pt.zl .* Δ.xl) ./ pt.xl) .* dat.lflag
    @timeit mpc.timer "Δzu" @. Δ.zu = ((ξxzu - pt.zu .* Δ.xu) ./ pt.xu) .* dat.uflag

    # IV. Set Δτ, Δκ to zero
    Δ.τ = zero(T)
    Δ.κ = zero(T)

    # Check Newton residuals
    # @printf "Newton residuals:\n"
    # @printf "|rp|   = %16.8e\n" norm(dat.A * Δ.x + mpc.regD .* Δ.y - ξp, Inf)
    # @printf "|rl|   = %16.8e\n" norm((Δ.x - Δ.xl) .* dat.lflag - ξl, Inf)
    # @printf "|ru|   = %16.8e\n" norm((Δ.x + Δ.xu) .* dat.uflag - ξu, Inf)
    # @printf "|rd|   = %16.8e\n" norm(-mpc.regP .* Δ.x + dat.A'Δ.y + Δ.zl - Δ.zu - ξd, Inf)
    # @printf "|rxzl| = %16.8e\n" norm(pt.zl .* Δ.xl + pt.xl .* Δ.zl - ξxzl, Inf)
    # @printf "|rxzu| = %16.8e\n" norm(pt.zu .* Δ.xu + pt.xu .* Δ.zu - ξxzu, Inf)

    return nothing
end

#=
"""
    max_step_length(x, dx)

Compute the maximum value `a ≥ 0` such that `x + a*dx ≥ 0`, where `x ≥ 0`.
"""
function max_step_length(x::Vector{T}, dx::Vector{T}) where{T}
    n = size(x, 1)
    n == size(dx, 1) || throw(DimensionMismatch())
    a = T(Inf)

    #=@inbounds=# for i in Base.OneTo(n)
        if dx[i] < zero(T)
            if (-x[i] / dx[i]) < a
                a = (-x[i] / dx[i])
            end
        end
    end
    return a
end

"""
    max_step_length(pt, δ)

Compute maximum length of homogeneous step.
"""
function max_step_length(pt::Point{T, Tv}, δ::Point{T, Tv}) where{T, Tv<:AbstractVector{T}}
    axl = max_step_length(pt.xl, δ.xl)
    axu = max_step_length(pt.xu, δ.xu)
    azl = max_step_length(pt.zl, δ.zl)
    azu = max_step_length(pt.zu, δ.zu)

    at = δ.τ < zero(T) ? (-pt.τ / δ.τ) : oneunit(T)
    ak = δ.κ < zero(T) ? (-pt.κ / δ.κ) : oneunit(T)
    
    α = min(one(T), axl, axu, azl, azu, at, ak)

    return α
end
=#


"""
    compute_higher_corrector!(Δc, mpc, γ, hx, hy, h0, Δ, α, β)

Compute higher-order corrected direction.

Requires the solution of one Newton system.

# Arguments
- `Δc`: Corrected search direction, modified in-place
- `mpc`: The MPC optimizer
- `γ`: 
- `hx, hy, h0`: Terms obtained from the preliminary augmented system solve
- `Δ`: Current predictor direction
- `α`: Maximum step length in predictor direction
- `β`: Relative threshold for centrality outliers
"""
function compute_higher_corrector!(Δc::Point{T, Tv},
    mpc::MPC{T, Tv}, γ::T,
    hx::Tv, hy::Tv, h0::T,
    Δ::Point{T, Tv}, α::T, β::T,
) where{T, Tv<:AbstractVector{T}}
    error("Higher corrector not implemented for MPC")
    # TODO: Sanity checks
    pt = mpc.pt
    dat = mpc.dat

    # Tentative step length
    α_ = min(one(T), T(2)*α)

    # Tentative cross products
    vl = ((pt.xl .+ α_ .* Δ.xl) .* (pt.zl .+ α_ .* Δ.zl)) .* dat.lflag
    vu = ((pt.xu .+ α_ .* Δ.xu) .* (pt.zu .+ α_ .* Δ.zu)) .* dat.uflag
    vt = (pt.τ   + α_  * Δ.τ)   * (pt.κ   + α_  * Δ.κ)

    # Compute target cross-products
    mu_l = β * pt.μ * γ
    mu_u = γ * pt.μ / β
    for i in 1:pt.n
        dat.lflag[i] || continue
        if vl[i] < mu_l
            vl[i] = mu_l - vl[i]
        elseif vl[i] > mu_u
            vl[i] = mu_u - vl[i]
        else
            vl[i] = zero(T)
        end
    end
    for i in 1:pt.n
        dat.uflag[i] || continue
        if vu[i] < mu_l
            vu[i] = mu_l - vu[i]
        elseif vu[i] > mu_u
            vu[i] = mu_u - vu[i]
        else
            vu[i] = zero(T)
        end
    end
    if vt < mu_l
        vt = mu_l - vt
    elseif vt > mu_u
        vt = mu_u - vt
    else
        vt = zero(T)
    end

    # Shift target cross-product to satisfy `v' * e = 0`
    δ = (sum(vl) + sum(vu) + vt) / (pt.p + 1)
    vl .-= δ
    vu .-= δ
    vt  -= δ

    # Compute corrector
    @timeit mpc.timer "Newton" solve_newton_system!(Δc, mpc, hx, hy, h0,
        # Right-hand sides
        tzeros(Tv, pt.m), tzeros(Tv, pt.n), tzeros(Tv, pt.n), tzeros(Tv, pt.n), zero(T),
        vl,
        vu,
        vt
    )

    # Update corrector
    Δc.x  .+= Δ.x
    Δc.xl .+= Δ.xl
    Δc.xu .+= Δ.xu
    Δc.y  .+= Δ.y
    Δc.zl .+= Δ.zl
    Δc.zu .+= Δ.zu
    Δc.τ   += Δ.τ
    Δc.κ   += Δ.κ

    # Compute corrected step-length
    αc = max_step_length(pt, Δc)
    return αc
end