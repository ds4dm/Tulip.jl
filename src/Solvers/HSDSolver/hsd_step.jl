"""
    compute_step!()

Compute next IP iterate for the HSD formulation.

# Arguments
- `model`: The optimization model
- `env`: Optimization environment
- `A`: Constraint matrix
- `F`: Factorization of the normal equations' matrix
- `b`: Right-hand side of primal constraints
- `c`: Primal linear objective term
- `uind, uval`: Indices and values of variables' upperbounds
- `θ, θwz`: Diagonal scaling terms
- `rp, ru, rd, rg`: Primal, dual and optimality residuals
- `x, w, y, s, z, t, k`: Primal-dual iterate
"""
function compute_step!(hsd::HSDSolver{Tv}, env::Env) where{Tv<:Real}

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
    θ .\= oneunit(Tv)

    # Update factorization
    update_linear_solver(hsd.ls, θ)

    # Search directions
    # Predictor
    Δ = Point(
        ncon, nvar, nupb,
        zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zeros(Tv, ncon), zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zero(Tv)
    )
    Δc = Point(
        ncon, nvar, nupb,
        zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zeros(Tv, ncon), zeros(Tv, nvar), zeros(Tv, nupb), zero(Tv),
        zero(Tv)
    )

    # Compute p, q, r, ρ from augmented system
    p = zeros(Tv, nvar)
    q = zeros(Tv, ncon)
    r = zeros(Tv, nupb)
    solve_augsys_hsd!(
        A, hsd.ls, θ, θwz, uind,
        p, q, r,
        b, c, uval
    )

    ρ = (pt.k / pt.t) - dot(c, p) + dot(b, q) - dot(uval, r)

    # Affine-scaling direction
    solve_newton_hsd!(
        A, hsd.ls, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        pt,
        Δ,
        # Right-hand side of Newton system
        res.rp, res.ru, res.rd, res.rg,
        -pt.x .* pt.s,
        -pt.w .* pt.z,
        -pt.t  * pt.k
    )

    # Step length for affine-scaling direction
    α = max_step_length(pt, Δ)
    γ = Tv((oneunit(Tv) - α)^2 * min(oneunit(Tv) - α, env.beta1))
    η = oneunit(Tv) - γ
    
    # Mehrothra corrector
    solve_newton_hsd!(
        A, hsd.ls, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        pt,
        Δ,
        # Right-hand side of Newton system
        η .* res.rp, η .* res.ru, η .* res.rd, η * res.rg,
        -pt.x .* pt.s .+ γ * pt.μ .- Δ.x .* Δ.s,
        -pt.w .* pt.z .+ γ * pt.μ .- Δ.w .* Δ.z,
        -pt.t  * pt.k  + γ * pt.μ  - Δ.t  * Δ.k
    )
    α = max_step_length(pt, Δ)

    # Extra corrections
    ncor = 0
    while (
        ncor < env.barrier_max_num_cor
        && α < Tv(0.999)    
    )
        α_ = α
        ncor += 1

        # TODO: Compute extra-corrector
        αc = compute_higher_corrector_hsd_!(
            ncon, nvar, nupb,
            η, γ,
            A, hsd.ls, b, c, uval, uind, θ, θwz,
            p, q, r, ρ,
            pt,
            Δ, α_,
            Δc,
            env.beta1, env.beta2, env.beta3, env.beta4,
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
    α *= Tv(0.99995)
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
    solve_newton_hsd!(...)

Solve Newton system
     A*dx           -b*dt                                   = ξp
     U*dx   +  dw   -u*dt                                   = ξu
                    -c*dt   +A'dy   +  ds   +U'dz           = ξd
    -c'dx                   +b'dy           -u'dz   -  dk   = ξg
     S*dx                           +X*ds                   = ξxs
            +Z*dw                          +W*dz            = ξwz
                     k*dt                           +t*dk   = ξtk

This version assumes that the augmented system below has been solved first.

# Arguments
- A: Matrix
- `F`: Factorization of the normal equations matrix
- `b`: Right-hand side of primal constraints
- `c`: Primal objective term
- `uind`: Indices of upper-bounded variables
- `uval`: Upper bound values
- `θ`: Diagonal scaling term
- `θwz`: Diagonal scaling term
- `p, q, r, ρ`: Terms obtained in the preliminary augmented system solve
- `x, w, y, s, z, t, k`: Primal and dual iterates
- `dx, dw, dy, ds, dz, dt, dk`: Primal and dual directions, modified in-place
- `ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk`: 
"""
function solve_newton_hsd!(
    A, ls, b, c, uind::Vector{Int}, uval, θ, θwz,
    p, q, r, ρ,
    pt::Point{Tv},
    Δ::Point{Tv},
    ξp, ξu, ξd, ξg::Tv, ξxs, ξwz, ξtk::Tv
) where{Tv<:Real}
    # Solve reduced newton system
    solve_augsys_hsd!(
        A, ls, θ, θwz, uind,
        Δ.x, Δ.y, Δ.z,
        ξp, 
        ξd .- (ξxs ./ pt.x),
        ξu .- (ξwz ./ pt.z)
    )

    # Compute Δτ
    Δ.t = (
        ξg + (ξtk / pt.t) + dot(c, Δ.x) - dot(b, Δ.y) + dot(uval, Δ.z)
    ) / ρ
    
    # Compute Δx, Δy, Δz
    Δ.x .+= Δ.t .* p
    Δ.y .+= Δ.t .* q
    Δ.z .+= Δ.t .* r
    # TODO: use functions below
    # axpy!(Δ.t, p, Δ.x)
    # axpy!(Δ.t, q, Δ.y)
    # axpy!(Δ.t, r, Δ.z)

    # Compute Δs, Δκ, Δw
    # ds  .= (rxs .- s  .* dx)   ./ x
    # dw  .= (rwz .- w  .* dz)   ./ z
    Δ.s  .= ξxs
    Δ.s .-= pt.s .* Δ.x
    Δ.s ./= pt.x

    Δ.w  .= ξwz
    Δ.w .-= pt.w .* Δ.z
    Δ.w ./= pt.z

    Δ.k = (ξtk - pt.k * Δ.t)  / pt.t

    return nothing
end


"""
    solve_augsys_hsd!(...)

Solve the augmented system below, and overwrite `dx, dy, dz` with the result.
    -(S/X)*dx  + A'dy  -    U'dz    = ξd
         A*dx                       = ξp
         U*dx          -(W/Z)*dz    = ξu

# Arguments
- `A`: Matrix
- `ls`: Linear solver for the augmented system
- `θ`: Diagonal scaling `1 ./ (s ./ x + U' (z ./ w) U)`
- `θwz`: Diagonal scaling `z ./ w`
- `uind`: Vector of indices of upper-bounded variables
- `dx, dy, dz`: Vectors of unknowns, modified in-place
- `ξp, ξd, ξu`: Right-hand-side vectors
"""
function solve_augsys_hsd!(
    A::AbstractMatrix{Tv}, ls::TLPLinearSolver{Tv},
    θ::Vector{Tv}, θwz::Vector{Tv}, uind::Vector{Int},
    dx::Vector{Tv}, dy::Vector{Tv}, dz::Vector{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}, ξu::Vector{Tv}
) where{Tv<:Real}
    m, n = size(A)
    
    # Set-up right-hand side
    ξd_ = copy(ξd)
    @views ξd_[uind] .-= (ξu .* θwz)

    solve_augmented_system!(dx, dy, ls, A, θ, ξp, ξd_)

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
    compute_higher_corrector(...)

Compute corrected Newton direction.

# Arguments
- `numvar, numvarub, numcon`: Number of variables and constraints
- `A, F, b, c, uval, uind`: Problem data
- `θ, θwz`: Diagonal scaling
- `p, q, r, ρ`: Vectors from initial augmented system
- `rp, ru, rd, rg`: Current primal and dual residuals
- `μ`: Current centering parameter
- `x, w, y, s, z, t, k`: Primal-dual iterate
- `dx, dw, dy, ds, dz, dt, dk`: Predictor direction (modified) 
- `α`: Maximum step length in predictor direction 
- `dxc, dwc, dyc, dsc, dzc, dtc, dkc`: Corrector direction (modified)
- `beta1, beta2, beta3, beta4`: Numerical parameters
"""
function compute_higher_corrector_hsd_!(
    ncon::Int, nvar::Int, nupb::Int,
    η::Tv, γ::Tv,
    A, ls, b, c, uval, uind::Vector{Int}, θ, θwz,
    p, q, r, ρ,
    pt::Point{Tv},
    Δ::Point{Tv}, α::Tv,
    Δc::Point{Tv},
    beta1::Float64, beta2::Float64, beta3::Float64, beta4::Float64,
) where{Tv<:Real}
    # TODO: Sanity checks

    # Tentative step length
    α_ = min(oneunit(Tv), Tv(2)*α)

    # Tentative cross products
    vx = (pt.x .+ α_ .* Δ.x) .* (pt.s .+ α_ .* Δ.s)
    vw = (pt.w .+ α_ .* Δ.w) .* (pt.z .+ α_ .* Δ.z)
    vt = (pt.t  + α_  * Δ.t)  * (pt.k  + α_  * Δ.k)

    # Compute target cross-products
    mu_l = Tv(beta4 * pt.μ * γ)
    mu_u = Tv(γ * pt.μ / beta4)
    #=@inbounds=# for i in 1:nvar
        if vx[i] < mu_l
            vx[i] = mu_l - vx[i]
        elseif vx[i] > mu_u
            vx[i] = mu_u - vx[i]
        else
            vx[i] = zero(Tv)
        end
    end
    #=@inbounds=# for i in 1:nupb
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
    δ = (sum(vx) + sum(vw) + vt) / (nvar + nupb + 1)
    vx .-= δ
    vw .-= δ
    vt  -= δ

    # Compute corrector
    solve_newton_hsd!(
        A, ls, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        pt,
        Δc,
        zeros(Tv, ncon), zeros(Tv, nupb), zeros(Tv, nvar), zero(Tv),
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