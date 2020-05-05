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
function compute_step!(hsd::HSDSolver{Tv}, params::Parameters{Tv}) where{Tv<:Real}

    # Names
    pt = hsd.pt
    res = hsd.res

    ncon = hsd.pt.m
    nvar = hsd.pt.n
    nupb = hsd.pt.p

    A = hsd.A
    b = hsd.b
    c = hsd.c
    uind = hsd.uind
    uval = hsd.uval

    # Compute scaling
    θ   = pt.s ./ pt.x
    θwz = pt.z ./ pt.w
    @views θ[uind] .+= θwz

    # Update regularizations
    rd_min = sqrt(eps(Tv))
    rp_min = sqrt(eps(Tv))
    rg_min = sqrt(eps(Tv))
    hsd.regP .= max.(rp_min, hsd.regP ./ 10)
    hsd.regD .= max.(rd_min, hsd.regD ./ 10)
    hsd.regG  = max( rg_min, hsd.regG  / 10)

    # Update factorization
    nbump = 0
    while nbump <= 3
        try
            KKT.update!(hsd.kkt, θ, hsd.regP, hsd.regD)
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
    nbump < 3 || throw(PosDefException(0))  # factorization could not be saved

    # Search directions
    # Predictor
    Δ  = Point{Tv}(ncon, nvar, nupb)
    Δc = Point{Tv}(ncon, nvar, nupb)

    # Compute hx, hy, hz from first augmented system solve
    hx = zeros(Tv, nvar)
    hy = zeros(Tv, ncon)
    hz = zeros(Tv, nupb)
    solve_augsys_hsd!(
        hx, hy, hz,
        hsd.kkt, θwz, uind,
        b, c, uval
    )
    # Recover h0 = ρg + κ / τ - c'hx + b'hy - u'hz
    h0 = (
        hsd.regG + pt.k / pt.t
        - dot(c, hx)
        + dot(b, hy)
        - dot(uval, hz)
    )

    # Affine-scaling direction
    solve_newton_hsd!(
        Δ,
        hsd.kkt, θwz, b, c, uind, uval,
        hx, hy, hz, h0, pt,
        # Right-hand side of Newton system
        res.rp, res.ru, res.rd, res.rg,
        -pt.x .* pt.s,
        -pt.w .* pt.z,
        -pt.t  * pt.k
    )

    # Step length for affine-scaling direction
    α = max_step_length(pt, Δ)
    γ = (oneunit(Tv) - α)^2 * min(oneunit(Tv) - α, params.BarrierGammaMin)
    η = oneunit(Tv) - γ
    
    # Mehrotra corrector
    solve_newton_hsd!(
        Δ,
        hsd.kkt, θwz, b, c, uind, uval,
        hx, hy, hz, h0, pt,
        # Right-hand side of Newton system
        η .* res.rp, η .* res.ru, η .* res.rd, η * res.rg,
        -pt.x .* pt.s .+ γ * pt.μ .- Δ.x .* Δ.s,
        -pt.w .* pt.z .+ γ * pt.μ .- Δ.w .* Δ.z,
        -pt.t  * pt.k  + γ * pt.μ  - Δ.t  * Δ.k
    )
    α = max_step_length(pt, Δ)

    # Extra corrections
    ncor = 0
    while ncor < params.BarrierCorrectionLimit && α < Tv(999 // 1000)
        α_ = α
        ncor += 1

        # TODO: Compute extra-corrector
        αc = compute_higher_corrector_hsd!(
            Δc, γ,
            hsd.kkt, θwz, b, c, uind, uval,
            hx, hy, hz, h0, pt,
            Δ, α_,
            params.BarrierCentralityOutlierThreshold,
        )
        if αc > α_
            # Use corrector
            Δ.x .= Δc.x
            Δ.w .= Δc.w
            Δ.y .= Δc.y
            Δ.s .= Δc.s
            Δ.z .= Δc.z
            Δ.t  = Δc.t
            Δ.k  = Δc.k
            α = αc
        end
        
        if αc < Tv(1.1) * α_
            break
        end

        # if (1.0 - η * αc) >= 0.9*(1.0 - η * α_)
        #     # not enough improvement, step correcting
        #     break
        # end
        
    end
    # Update current iterate
    α *= params.BarrierStepDampFactor
    pt.x .+= α .* Δ.x
    pt.w .+= α .* Δ.w
    pt.t  += α  * Δ.t
    pt.y .+= α .* Δ.y
    pt.s .+= α .* Δ.s
    pt.z .+= α .* Δ.z
    pt.k  += α  * Δ.k
    update_mu!(pt)
    
    return nothing
end


"""
    solve_newton_hsd!

Solve the Newton system for `dx, dw, dy, ds, dz, dτ, dκ`
```
-Rp*dx        + A'dy  +   ds - U'dz -  c*dτ        = ξd
  A*dx        + Rd*dy               -  b*dτ        = ξp
  U*dx +   dw                       -  u*dτ        = ξu
 -c'dx        + b'dy         - u'dz + Rg*dτ -   dκ = ξg
  S*dx                + X*ds                       = ξxs
         Z*dw                + W*dz                = ξwz
                                       k*dτ + τ*dκ = ξτκ
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
function solve_newton_hsd!(
    Δ::Point{Tv},
    kkt, θwz, b, c, uind::Vector{Int}, uval,
    hx::Vector{Tv}, hy::Vector{Tv}, hz::Vector{Tv}, h0::Tv, pt::Point{Tv},
    ξp::Vector{Tv}, ξu::Vector{Tv}, ξd::Vector{Tv}, ξg::Tv,
    ξxs::Vector{Tv}, ξwz::Vector{Tv}, ξtk::Tv
) where{Tv<:Real}

    # Solve reduced newton system
    solve_augsys_hsd!(
        Δ.x, Δ.y, Δ.z,
        kkt, θwz, uind,
        ξp, ξd .- (ξxs ./ pt.x), ξu .- (ξwz ./ pt.z)
    )

    # Compute Δτ
    Δ.t = (
        ξg + (ξtk / pt.t)
        + dot(c, Δ.x)
        - dot(b, Δ.y)
        + dot(uval, Δ.z)
    ) / h0
    
    # Compute Δx, Δy, Δz
    Δ.x .+= Δ.t .* hx
    Δ.y .+= Δ.t .* hy
    Δ.z .+= Δ.t .* hz
    # TODO: use functions below
    # axpy!(Δ.t, hx, Δ.x)
    # axpy!(Δ.t, hy, Δ.y)
    # axpy!(Δ.t, hz, Δ.z)

    # Recover Δs, Δw and Δk
    # Δs = (ξxs - s .* Δx) ./ x
    Δ.s  .= ξxs
    Δ.s .-= pt.s .* Δ.x
    Δ.s ./= pt.x

    # Δw = (ξwz - w .* Δx) ./ z
    Δ.w  .= ξwz
    Δ.w .-= pt.w .* Δ.z
    Δ.w ./= pt.z

    Δ.k = (ξtk - pt.k * Δ.t)  / pt.t

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
function max_step_length(x::Vector{Tv}, dx::Vector{Tv}) where{Tv<:Real}
    n = size(x, 1)
    n == size(dx, 1) || throw(DimensionMismatch())
    a = Tv(Inf)

    #=@inbounds=# for i in Base.OneTo(n)
        if dx[i] < zero(Tv)
            if (-x[i] / dx[i]) < a
                a = (-x[i] / dx[i])
            end
        end
    end
    return a
end

function max_step_length(pt::Point{Tv}, δ::Point{Tv}) where{Tv<:Real}
    ax = max_step_length(pt.x, δ.x)
    aw = max_step_length(pt.w, δ.w)
    as = max_step_length(pt.s, δ.s)
    az = max_step_length(pt.z, δ.z)

    at = δ.t < zero(Tv) ? (-pt.t / δ.t) : oneunit(Tv)
    ak = δ.k < zero(Tv) ? (-pt.k / δ.k) : oneunit(Tv)
    
    α = min(oneunit(Tv), ax, aw, as, az, at, ak)

    return α
end


"""
    compute_higher_corrector_hsd!

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
function compute_higher_corrector_hsd!(
    Δc::Point{Tv},
    γ::Tv,
    kkt, θwz, b, c, uind::Vector{Int}, uval,
    hx, hy, hz, h0, pt::Point{Tv},
    Δ::Point{Tv}, α::Tv,
    β::Tv,
) where{Tv<:Real}
    # TODO: Sanity checks
    
    # Tentative step length
    α_ = min(oneunit(Tv), Tv(2)*α)

    # Tentative cross products
    vx = (pt.x .+ α_ .* Δ.x) .* (pt.s .+ α_ .* Δ.s)
    vw = (pt.w .+ α_ .* Δ.w) .* (pt.z .+ α_ .* Δ.z)
    vt = (pt.t  + α_  * Δ.t)  * (pt.k  + α_  * Δ.k)

    # Compute target cross-products
    mu_l = β * pt.μ * γ
    mu_u = γ * pt.μ / β
    for i in 1:pt.n
        if vx[i] < mu_l
            vx[i] = mu_l - vx[i]
        elseif vx[i] > mu_u
            vx[i] = mu_u - vx[i]
        else
            vx[i] = zero(Tv)
        end
    end
    for i in 1:pt.p
        if vw[i] < mu_l
            vw[i] = mu_l - vw[i]
        elseif vx[i] > mu_u
            vw[i] = mu_u - vw[i]
        else
            vw[i] = zero(Tv)
        end
    end
    if vt < mu_l
        vt = mu_l - vt
    elseif vt > mu_u
        vt = mu_u - vt
    else
        vt = zero(Tv)
    end

    # Shift target cross-product to satisfy `v' * e = 0`
    δ = (sum(vx) + sum(vw) + vt) / (pt.n + pt.p + 1)
    vx .-= δ
    vw .-= δ
    vt  -= δ

    # Compute corrector
    solve_newton_hsd!(
        Δc,
        kkt, θwz, b, c, uind, uval,
        hx, hy, hz, h0, pt,
        # Right-hand sides
        zeros(Tv, pt.m), zeros(Tv, pt.p), zeros(Tv, pt.n), zero(Tv),
        vx,
        vw,
        vt
    )

    # Update corrector
    Δc.x .+= Δ.x
    Δc.w .+= Δ.w
    Δc.y .+= Δ.y
    Δc.s .+= Δ.s
    Δc.z .+= Δ.z
    Δc.t  += Δ.t
    Δc.k  += Δ.k

    # Compute corrected step-length
    αc = max_step_length(pt, Δc)
    return αc

end