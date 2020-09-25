"""
    compute_step!()

Compute next IP iterate for the HSD formulation.

# Arguments
- `model`: The optimization model
- `params`: Optimization parameters
- `A`: Constraint matrix
- `F`: Factorization of the normal equations' matrix
- `b`: Right-hand side of primal constraints
- `c`: Primal linear objective term
- `uind, uval`: Indices and values of variables' upperbounds
- `θ, θwz`: Diagonal scaling terms
- `rp, ru, rd, rg`: Primal, dual and optimality residuals
- `x, w, y, s, z, t, k`: Primal-dual iterate
"""
function compute_step!(hsd::HSDSolver{T, Tv, Tk},
    dat::IPMData{T, Tv, Tb, Ta}, params::Parameters{T}
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}, Tk<:AbstractKKTSolver{T}}

    # Names
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

    # Check that no NaN values appeared
    @assert !any(isnan.(θl))
    @assert !any(isnan.(θu))
    @assert !any(isnan.(θinv))

    # Update regularizations
    hsd.regP .= max.(params.BarrierPRegMin, hsd.regP ./ 10)
    hsd.regD .= max.(params.BarrierDRegMin, hsd.regD ./ 10)
    hsd.regG  = max( params.BarrierPRegMin, hsd.regG  / 10)

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
    Δ  = Point{T, Tv}(m, n, p)
    Δc = Point{T, Tv}(m, n, p)

    # Compute hx, hy, hz from first augmented system solve
    hx = tzeros(Tv, n)
    hy = tzeros(Tv, m)
    ξ_ = dat.c - ((pt.zl ./ pt.xl) .* dat.l) .* dat.lflag - ((pt.zu ./ pt.xu) .* dat.u) .* dat.uflag
    KKT.solve!(hx, hy, hsd.kkt, dat.b, ξ_)

    # Recover h0 = ρg + κ / τ - c'hx + b'hy - u'hz
    h0 = (
        hsd.regG + pt.κ / pt.τ
        + (dot(dat.l .* dat.lflag, (dat.l .* θl) .* dat.lflag)
        + dot(dat.u .* dat.uflag, (dat.u .* θu) .* dat.uflag)
        - dot(c + (θl .* dat.l) .* dat.lflag + (θu .* dat.u) .* dat.uflag, hx))
        + dot(b, hy)
    )

    # Affine-scaling direction
    @timeit hsd.timer "Newton" solve_newton_system!(Δ, hsd, dat, hx, hy, h0,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd, res.rg,
        -(pt.xl .* pt.zl) .* dat.lflag,
        -(pt.xu .* pt.zu) .* dat.uflag,
        -pt.τ  * pt.κ
    )

    # Step length for affine-scaling direction
    α = max_step_length(pt, Δ)
    # @info "Affine step length: $α"
    γ = (one(T) - α)^2 * min(one(T) - α, params.BarrierGammaMin)
    η = one(T) - γ
    
    # Mehrotra corrector
    @timeit hsd.timer "Newton" solve_newton_system!(Δ, hsd, dat, hx, hy, h0,
        # Right-hand side of Newton system
        η .* res.rp, η .* res.rl, η .* res.ru, η .* res.rd, η * res.rg,
        (-pt.xl .* pt.zl .+ γ * pt.μ .- Δ.xl .* Δ.zl) .* dat.lflag,
        (-pt.xu .* pt.zu .+ γ * pt.μ .- Δ.xu .* Δ.zu) .* dat.uflag,
        -pt.τ  * pt.κ  + γ * pt.μ  - Δ.τ  * Δ.κ
    )
    α = max_step_length(pt, Δ)
    # @info "Corrector step length: $α"

    # Extra corrections
    ncor = 0
    while ncor < params.BarrierCorrectionLimit && α < T(999 // 1000)
        α_ = α
        ncor += 1

        # TODO: Compute extra-corrector
        αc = compute_higher_corrector_hsd!(Δc, 
            hsd, dat, γ, 
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
        
    end
    # Update current iterate
    α *= params.BarrierStepDampFactor
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
    solve_newton_system!(Δ, hsd::HSDSolver, dat::IPMData)

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
- `kkt`: Linear solver for the augmented system
- `θwz`: Diagonal scaling term
- `b`: Right-hand side of primal constraints
- `c`: Primal objective term
- `uind`: Indices of upper-bounded variables
- `uval`: Upper bound values
- `hx, hy, hz, h0`: Terms obtained in the preliminary augmented system solve
- `pt`: Current primal-dual iterate
- `ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk`: Right-hand side vectors
"""
function solve_newton_system!(Δ::Point{T, Tv},
    hsd::HSDSolver{T, Tv, Tk},
    dat::IPMData{T, Tv, Tb, Ta},
    # Information from initial augmented system solve
    hx::Tv, hy::Tv, h0::T,
    # Right-hand side
    ξp::Tv, ξl::Tv, ξu::Tv, ξd::Tv, ξg::T, ξxzl::Tv, ξxzu::Tv, ξtk::T
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}, Tk<:AbstractKKTSolver{T}}

    # Check that no NaN values are in the RHS
    # @assert !any(isnan.(ξp))
    # @assert !any(isnan.(ξl))
    # @assert !any(isnan.(ξu))
    # @assert !any(isnan.(ξd))
    # @assert !any(isnan.(ξg))
    # @assert !any(isnan.(ξxzl))
    # @assert !any(isnan.(ξxzu))
    # @assert !any(isnan.(ξtk))

    pt = hsd.pt

    # I. Solve augmented system
    @timeit hsd.timer "ξd_"  begin
        ξd_ = copy(ξd)
        @. ξd_ += -((ξxzl + pt.zl .* ξl) ./ pt.xl) .* dat.lflag + ((ξxzu - pt.zu .* ξu) ./ pt.xu) .* dat.uflag
    end
    @timeit hsd.timer "KKT" KKT.solve!(Δ.x, Δ.y, hsd.kkt, ξp, ξd_)

    # @assert !any(isnan.(Δ.x))
    # @assert !any(isnan.(Δ.y))

    # II. Recover Δτ, Δx, Δy
    # Compute Δτ
    @timeit hsd.timer "ξg_" ξg_ = (ξg + ξtk / pt.τ 
        - dot((ξxzl ./ pt.xl) .* dat.lflag, dat.l .* dat.lflag)            # l'(Xl)^-1 * ξxzl
        + dot((ξxzu ./ pt.xu) .* dat.uflag, dat.u .* dat.uflag) 
        - dot(((pt.zl ./ pt.xl) .* ξl) .* dat.lflag, dat.l .* dat.lflag) 
        - dot(((pt.zu ./ pt.xu) .* ξu) .* dat.uflag, dat.u .* dat.uflag)  # 
    )

    # @assert !isnan(ξg_)

    @timeit hsd.timer "Δτ" Δ.τ = (
        ξg_
        + dot(dat.c
            + ((pt.zl ./ pt.xl) .* dat.l) .* dat.lflag
            + ((pt.zu ./ pt.xu) .* dat.u) .* dat.uflag
        , Δ.x)
        - dot(dat.b, Δ.y)
    ) / h0
    
    # @assert !isnan(Δ.τ)

    # Compute Δx, Δy
    @timeit hsd.timer "Δx" Δ.x .+= Δ.τ .* hx
    @timeit hsd.timer "Δy" Δ.y .+= Δ.τ .* hy

    # @assert !any(isnan.(Δ.x))
    # @assert !any(isnan.(Δ.y))

    # III. Recover Δxl, Δxu
    @timeit hsd.timer "Δxl" begin
        @. Δ.xl = -ξl + Δ.x - Δ.τ .* (dat.l .* dat.lflag)
        Δ.xl .*= dat.lflag
    end
    @timeit hsd.timer "Δxu" begin
        @. Δ.xu =  ξu - Δ.x + Δ.τ .* (dat.u .* dat.uflag)
        Δ.xu .*= dat.uflag
    end

    # @assert !any(isnan.(Δ.xl))
    # @assert !any(isnan.(Δ.xu))

    # IV. Recover Δzl, Δzu
    @timeit hsd.timer "Δzl" @. Δ.zl = ((ξxzl - pt.zl .* Δ.xl) ./ pt.xl) .* dat.lflag
    @timeit hsd.timer "Δzu" @. Δ.zu = ((ξxzu - pt.zu .* Δ.xu) ./ pt.xu) .* dat.uflag

    # V. Recover Δκ
    Δ.κ = (ξtk - pt.κ * Δ.τ) / pt.τ
    # @assert !any(isnan.(Δ.κ))

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

    # @assert all(iszero.(Δ.xl .* .!(dat.lflag)))
    # @assert all(iszero.(Δ.xu .* .!(dat.uflag)))
    # @assert all(iszero.(Δ.zl .* .!(dat.lflag)))
    # @assert all(iszero.(Δ.zu .* .!(dat.uflag)))

    return nothing
end


"""
    solve_augsys_hsd!

Solve the augmented system for `dx, dy, dz`
```
-(S/X + Rp)*dx +  A'dy     - U'dz = ξd
          A*dx + Rd*dy            = ξp
          U*dx          -(Z/W)*dz = ξu
```

# Arguments
- `dx, dy, dz`: Vectors of unknowns, modified in-place
- `kkt`: Linear solver for the augmented system
- `θwz`: Diagonal scaling `z ./ w`
- `uind`: Vector of indices of upper-bounded variables
- `ξp, ξd, ξu`: Right-hand-side vectors
"""
function solve_augsys_hsd!(
    dx::Vector{Tv}, dy::Vector{Tv}, dz::Vector{Tv},
    kkt::AbstractKKTSolver{Tv},
    θwz::Vector{Tv}, uind::Vector{Int},
    ξp::Vector{Tv}, ξd::Vector{Tv}, ξu::Vector{Tv}
) where{Tv<:Real}
    
    # Set-up right-hand side
    ξd_ = copy(ξd)
    @views ξd_[uind] .-= (ξu .* θwz)

    KKT.solve!(dx, dy, kkt, ξp, ξd_)

    # Recover dz
    dz .= zero(Tv)
    dz .-= ξu
    @views dz .+= dx[uind]
    dz .*= θwz

    return nothing
end


"""
    max_step_length(x, dx)

Compute maximum step size `a > 0` such that `x + a*dx >= 0`, where `x` is a 
    non-negative vector.

    max_step_length(...)

Compute maximum length of homogeneous step.
"""
function max_step_length(x::Vector{T}, dx::Vector{T}) where{T<:Real}
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
    compute_higher_corrector!(Δc, hsd, dat, γ, hx, hy, h0, Δ, α, β)

Compute higher-order corrected direction.

Requires the solution of one Newton system.

# Arguments
- `Δc`: Corrected search direction, modified in-place
- `γ`
- `kkt`: Linear solver for the augmented system
- `θwz`: Diagonal scaling term
- `b`: Right-hand side of primal constraints
- `c`: Primal objective term
- `uind`: Indices of upper-bounded variables
- `uval`: Upper bound values
- `hx, hy, hz, h0`: Terms obtained in the preliminary augmented system solve
- `regP, regD`: Primal and dual regularizers
- `x̄, ȳ`: Primal and dual proximal points
- `pt`: Current primal-dual iterate
- `Δ`: Current predictor direction
- `α`: Maximum step length in predictor direction 
- `β`: Relative threshold for centrality outliers
"""
function compute_higher_corrector_hsd!(Δc::Point{T, Tv},
    hsd::HSDSolver{T, Tv, Tk}, dat::IPMData{T, Tv, Tb, Ta}, γ::T,
    hx::Tv, hy::Tv, h0::T,
    Δ::Point{T, Tv}, α::T, β::T,
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}, Tk<:AbstractKKTSolver{T}}
    # TODO: Sanity checks
    pt = hsd.pt

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
    @timeit hsd.timer "Newton" solve_newton_system!(Δc, hsd, dat, hx, hy, h0,
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