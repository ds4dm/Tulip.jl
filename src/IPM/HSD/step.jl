"""
    compute_step!(hsd, params)

Compute next IP iterate for the HSD formulation.

# Arguments
- `hsd`: The HSD optimizer model
- `params`: Optimization parameters
"""
function compute_step!(hsd::HSD{T, Tv}, params::IPMOptions{T}) where{T, Tv<:AbstractVector{T}}

    # Names
    dat = hsd.dat
    pt = hsd.pt
    res = hsd.res

    m, n, p = pt.m, pt.n, pt.p

    A = dat.A
    b = dat.b
    c = dat.c

    # Compute scaling
    θl = (pt.zl ./ pt.xl) .* dat.lflag
    θu = (pt.zu ./ pt.xu) .* dat.uflag
    θinv = θl .+ θu

    # Update regularizations
    hsd.regP .= max.(params.PRegMin, hsd.regP ./ 10)
    hsd.regD .= max.(params.DRegMin, hsd.regD ./ 10)
    hsd.regG  = max( params.PRegMin, hsd.regG  / 10)

    # Update factorization
    nbump = 0
    while nbump <= 3
        try
            @timeit hsd.timer "Factorization" KKT.update!(hsd.kkt, θinv, hsd.regP, hsd.regD)
            break
        catch err
            isa(err, PosDefException) || isa(err, ZeroPivotException) || rethrow(err)

            # Increase regularization
            hsd.regD .*= 100
            hsd.regP .*= 100
            hsd.regG  *= 100
            nbump += 1
            @warn "Increase regularizations to $(hsd.regG)"
        end
    end
    # TODO: throw a custom error for numerical issues
    nbump < 3 || throw(PosDefException(0))  # factorization could not be saved

    # Search directions
    # Predictor
    Δ  = Point{T, Tv}(m, n, p, hflag=true)
    Δc = Point{T, Tv}(m, n, p, hflag=true)

    # Compute hx, hy, hz from first augmented system solve
    hx = tzeros(Tv, n)
    hy = tzeros(Tv, m)
    ξ_ = @. (dat.c - ((pt.zl / pt.xl) * dat.l) * dat.lflag - ((pt.zu / pt.xu) * dat.u) * dat.uflag)
    @timeit hsd.timer "Newton" begin
      @timeit hsd.timer "KKT" KKT.solve!(hx, hy, hsd.kkt, dat.b, ξ_)
    end

    dot_buf = buffer_for_dot_weighted_sum(T)

    # Recover h0 = ρg + κ / τ - c'hx + b'hy - u'hz
    # Some of the summands may take large values,
    # so care must be taken for numerical stability
    h0 = buffered_dot_weighted_sum!!(
        dot_buf,
        (
            (dat.l[dat.lflag], (dat.l .* θl)[dat.lflag]),
            (dat.u[dat.uflag], (dat.u .* θu)[dat.uflag]),
            ((@. (c + (θl * dat.l) * dat.lflag + (θu * dat.u) * dat.uflag)), hx),
            (b, hy),
        ),
        (
            1, 1, -1, 1,
        ),
    ) + pt.κ / pt.τ + hsd.regG

    # Affine-scaling direction
    @timeit hsd.timer "Newton" solve_newton_system!(Δ, hsd, hx, hy, h0,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd, res.rg,
        .-(pt.xl .* pt.zl) .* dat.lflag,
        .-(pt.xu .* pt.zu) .* dat.uflag,
        .-pt.τ  * pt.κ
    )

    # Step length for affine-scaling direction
    α = max_step_length(pt, Δ)
    γ = (one(T) - α)^2 * min(one(T) - α, params.GammaMin)
    η = one(T) - γ

    # Mehrotra corrector
    @timeit hsd.timer "Newton" solve_newton_system!(Δ, hsd, hx, hy, h0,
        # Right-hand side of Newton system
        η .* res.rp, η .* res.rl, η .* res.ru, η .* res.rd, η * res.rg,
        (.-pt.xl .* pt.zl .+ γ * pt.μ .- Δ.xl .* Δ.zl) .* dat.lflag,
        (.-pt.xu .* pt.zu .+ γ * pt.μ .- Δ.xu .* Δ.zu) .* dat.uflag,
        -pt.τ  * pt.κ  + γ * pt.μ  - Δ.τ  * Δ.κ
    )
    α = max_step_length(pt, Δ)

    # Extra corrections
    ncor = 0
    while ncor < params.CorrectionLimit && α < T(999 // 1000)
        α_ = α
        ncor += 1

        # Compute extra-corrector
        αc = compute_higher_corrector!(Δc,
            hsd, γ,
            hx, hy, h0,
            Δ, α_, params.CentralityOutlierThreshold
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

    end

    # Update current iterate
    α *= params.StepDampFactor
    pt.x  .+= α .* Δ.x
    pt.xl .+= α .* Δ.xl
    pt.xu .+= α .* Δ.xu
    pt.y  .+= α .* Δ.y
    pt.zl .+= α .* Δ.zl
    pt.zu .+= α .* Δ.zu
    pt.τ   += α  * Δ.τ
    pt.κ   += α  * Δ.κ
    update_mu!(pt)

    return nothing
end


"""
    solve_newton_system!(Δ, hsd, hx, hy, h0, ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk)

Solve the Newton system
```math
\\begin{bmatrix}
    A & & & R_{d} & & & -b\\\\
    I & -I & & & & & -l\\\\
    I & & -I & & & & -u\\\\
    -R_{p} & & & A^{T} & I & -I & -c\\\\
    -c^{T} & & & b^{T} & l^{T} & -u^{T} & ρ_{g} & -1\\\\
    & Z_{l} & & & X_{l}\\\\
    & & Z_{u} & & & X_{u}\\\\
    &&&&&& κ & τ
\\end{bmatrix}
\\begin{bmatrix}
    Δ x\\\\
    Δ x_{l}\\\\
    Δ x_{u}\\\\
    Δ y\\\\
    Δ z_{l} \\\\
    Δ z_{u} \\\\
    Δ τ\\\\
    Δ κ\\\\
\\end{bmatrix}
=
\\begin{bmatrix}
    ξ_p\\\\
    ξ_l\\\\
    ξ_u\\\\
    ξ_d\\\\
    ξ_g\\\\
    ξ_{xz}^{l}\\\\
    ξ_{xz}^{u}\\\\
    ξ_tk
\\end{bmatrix}
```

# Arguments
- `Δ`: Search direction, modified
- `hsd`: The HSD optimizer
- `hx, hy, hz, h0`: Terms obtained in the preliminary augmented system solve
- `ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk`: Right-hand side vectors
"""
function solve_newton_system!(Δ::Point{T, Tv},
    hsd::HSD{T, Tv},
    # Information from initial augmented system solve
    hx::Tv, hy::Tv, h0::T,
    # Right-hand side
    ξp::Tv, ξl::Tv, ξu::Tv, ξd::Tv, ξg::T, ξxzl::Tv, ξxzu::Tv, ξtk::T
) where{T, Tv<:AbstractVector{T}}

    pt = hsd.pt
    dat = hsd.dat

    # I. Solve augmented system
    @timeit hsd.timer "ξd_"  begin
        ξd_ = copy(ξd)
        @. ξd_ += -((ξxzl + pt.zl .* ξl) ./ pt.xl) .* dat.lflag + ((ξxzu - pt.zu .* ξu) ./ pt.xu) .* dat.uflag
    end
    @timeit hsd.timer "KKT" KKT.solve!(Δ.x, Δ.y, hsd.kkt, ξp, ξd_)

    dot_buf = buffer_for_dot_weighted_sum(T)

    # II. Recover Δτ, Δx, Δy
    # Compute Δτ
    @timeit hsd.timer "ξg_" ξg_ = ξg + ξtk / pt.τ +
        buffered_dot_weighted_sum!!(
            dot_buf,
            (
                ((ξxzl ./ pt.xl)[dat.lflag], dat.l[dat.lflag]),            # l'(Xl)^-1 * ξxzl
                ((ξxzu ./ pt.xu)[dat.uflag], dat.u[dat.uflag]),
                (((pt.zl ./ pt.xl) .* ξl)[dat.lflag], dat.l[dat.lflag]),
                (((pt.zu ./ pt.xu) .* ξu)[dat.uflag], dat.u[dat.uflag]),
            ),
            (
                -1, 1, -1, -1,
            ),
        )

    @timeit hsd.timer "Δτ" Δ.τ = (
        ξg_ +
        buffered_dot_weighted_sum!!(
            dot_buf,
            (
                (
                    (@. (
                        dat.c +
                        ((pt.zl / pt.xl) * dat.l) * dat.lflag +
                        ((pt.zu / pt.xu) * dat.u) * dat.uflag)),
                    Δ.x,
                ),
                (dat.b, Δ.y),
            ),
            (
                1, -1,
            ),
        )
    ) / h0


    # Compute Δx, Δy
    @timeit hsd.timer "Δx" Δ.x .+= Δ.τ .* hx
    @timeit hsd.timer "Δy" Δ.y .+= Δ.τ .* hy

    # III. Recover Δxl, Δxu
    @timeit hsd.timer "Δxl" begin
        @. Δ.xl = (-ξl + Δ.x - Δ.τ .* (dat.l .* dat.lflag)) * dat.lflag
    end
    @timeit hsd.timer "Δxu" begin
        @. Δ.xu = ( ξu - Δ.x + Δ.τ .* (dat.u .* dat.uflag)) * dat.uflag
    end

    # IV. Recover Δzl, Δzu
    @timeit hsd.timer "Δzl" @. Δ.zl = ((ξxzl - pt.zl .* Δ.xl) ./ pt.xl) .* dat.lflag
    @timeit hsd.timer "Δzu" @. Δ.zu = ((ξxzu - pt.zu .* Δ.xu) ./ pt.xu) .* dat.uflag

    # V. Recover Δκ
    Δ.κ = (ξtk - pt.κ * Δ.τ) / pt.τ

    # Check Newton residuals
    # @printf "Newton residuals:\n"
    # @printf "|rp|   = %16.8e\n" norm(dat.A * Δ.x + hsd.regD .* Δ.y - dat.b .* Δ.τ - ξp, Inf)
    # @printf "|rl|   = %16.8e\n" norm((Δ.x - Δ.xl - (dat.l .* Δ.τ)) .* dat.lflag - ξl, Inf)
    # @printf "|ru|   = %16.8e\n" norm((Δ.x + Δ.xu - (dat.u .* Δ.τ)) .* dat.uflag - ξu, Inf)
    # @printf "|rd|   = %16.8e\n" norm(-hsd.regP .* Δ.x + dat.A'Δ.y + Δ.zl - Δ.zu - dat.c .* Δ.τ - ξd, Inf)
    # @printf "|rg|   = %16.8e\n" norm(-dat.c'Δ.x + dat.b'Δ.y + dot(dat.l .* dat.lflag, Δ.zl) - dot(dat.u .* dat.uflag, Δ.zu) + hsd.regG * Δ.τ - Δ.κ - ξg, Inf)
    # @printf "|rxzl| = %16.8e\n" norm(pt.zl .* Δ.xl + pt.xl .* Δ.zl - ξxzl, Inf)
    # @printf "|rxzu| = %16.8e\n" norm(pt.zu .* Δ.xu + pt.xu .* Δ.zu - ξxzu, Inf)
    # @printf "|rtk|  = %16.8e\n" norm(pt.κ * Δ.τ + pt.τ * Δ.κ - ξtk, Inf)

    return nothing
end


"""
    max_step_length(x, dx)

Compute the maximum value `a ≥ 0` such that `x + a*dx ≥ 0`, where `x ≥ 0`.
"""
function max_step_length(x::Vector{T}, dx::Vector{T}) where{T}
    n = size(x, 1)
    n == size(dx, 1) || throw(DimensionMismatch())
    a = T(Inf)

    @inbounds for i in 1:n
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


"""
    compute_higher_corrector!(Δc, hsd, γ, hx, hy, h0, Δ, α, β)

Compute higher-order corrected direction.

Requires the solution of one Newton system.

# Arguments
- `Δc`: Corrected search direction, modified in-place
- `hsd`: The HSD optimizer
- `γ`:
- `hx, hy, h0`: Terms obtained from the preliminary augmented system solve
- `Δ`: Current predictor direction
- `α`: Maximum step length in predictor direction
- `β`: Relative threshold for centrality outliers
"""
function compute_higher_corrector!(Δc::Point{T, Tv},
    hsd::HSD{T, Tv}, γ::T,
    hx::Tv, hy::Tv, h0::T,
    Δ::Point{T, Tv}, α::T, β::T,
) where{T, Tv<:AbstractVector{T}}
    # TODO: Sanity checks
    pt = hsd.pt
    dat = hsd.dat

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
    @timeit hsd.timer "Newton" solve_newton_system!(Δc, hsd, hx, hy, h0,
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
