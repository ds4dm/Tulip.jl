"""
    make_step!()

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
function make_step!(
    ::Val{1},
    model::Model, env::TulipEnv,
    A, F, b, c, uind::Vector{Int}, uval::Vector{Float64},
    θ, θwz, μ,
    rp, ru, rd, rg::RefValue,
    x, w, y, s, z, t::RefValue, k::RefValue
)

    # Search directions
    # Predictor
    dx = zeros(model.num_var)
    dw = zeros(model.n_var_ub)
    dy = zeros(model.num_constr)
    ds = zeros(model.num_var)
    dz = zeros(model.n_var_ub)
    dt = Ref(0.0)
    dk = Ref(0.0)

    # Corrector
    dxc = zeros(model.num_var)
    dwc = zeros(model.n_var_ub)
    dyc = zeros(model.num_constr)
    dsc = zeros(model.num_var)
    dzc = zeros(model.n_var_ub)
    dtc = Ref(0.0)
    dkc = Ref(0.0)

    # Compute p, q, r, ρ from augmented system
    p = zeros(model.num_var)
    q = zeros(model.num_constr)
    r = zeros(model.n_var_ub)
    solve_augsys_hsd!(
        A, F, θ, θwz, uind,
        p, q, r,
        b, c, uval
    )
    ρ = (k.x / t.x) - dot(c, p) + dot(b, q) - dot(uval, r)

    # Affine-scaling direction
    solve_newton_hsd!(
        A, F, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        x, w, y, s, z, t, k,
        dx, dw, dy, ds, dz, dt, dk,
        # Right-hand side of Newton system
        rp, ru, rd, rg.x,
        -x  .* s,
        -w  .* z,
        -t.x * k.x
    )

    # Step length for affine-scaling direction
    α = max_step_length(x, w, y, s, z, t, k, dx, dw, dy, ds, dz, dt, dk)
    γ = (1-α)^2 * min(1-α, env.beta1.val)
    η = 1.0 - γ
    
    # Mehrothra corrector
    solve_newton_hsd!(
        A, F, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        x, w, y, s, z, t, k,
        dx, dw, dy, ds, dz, dt, dk,
        # Right-hand side of Newton system
        η .* rp, η .* ru, η .* rd, η * rg.x,
        -x  .* s   .+ γ * μ.x .- dx  .* ds,
        -w  .* z   .+ γ * μ.x .- dw  .* dz,
        -t.x * k.x  + γ * μ.x  - dt.x * dk.x
    )
    α = max_step_length(
        x, w, y, s, z, t, k,
        dx, dw, dy, ds, dz, dt, dk
    )

    # Extra corrections
    ncor = 0
    while (
        ncor < env.barrier_max_num_cor.val
        && α < 0.999    
    )
        α_ = α
        ncor += 1

        # TODO: Compute extra-corrector
        αc = compute_higher_corrector_hsd_!(
            model.num_var, model.n_var_ub, model.num_constr,
            η, γ,
            A, F, b, c, uval, uind, θ, θwz,
            p, q, r, ρ,
            x, w, y, s, z, t, k,
            rp, ru, rd, rg, μ,
            dx, dw, dy, ds, dz, dt, dk, α_,
            dxc, dwc, dyc, dsc, dzc, dtc, dkc,
            env.beta1.val, env.beta2.val, env.beta3.val, env.beta4.val,
        )
        if αc > α_
            # Use corrector
            dx   .= dxc
            dw   .= dwc
            dy   .= dyc
            ds   .= dsc
            dz   .= dzc
            dt.x  = dtc.x
            dk.x  = dkc.x
            α = αc
        end
        
        if αc < 1.1 * α_
            break
        end

        # if (1.0 - η * αc) >= 0.9*(1.0 - η * α_)
        #     # not enough improvement, step correcting
        #     break
        # end
        
    end
    # Update current iterate
    α *= 0.99995
    disable_sigint() do
        x  .+= α .* dx
        w  .+= α .* dw
        t.x += α  * dt.x
        y  .+= α .* dy
        s  .+= α .* ds
        z  .+= α .* dz
        k.x += α  * dk.x
    end
    
    return nothing
end

"""
    solve_newton_hsd!(...)

Solve Newton system
     A*dx           -b*dt                                   = rp
     U*dx   +  dw   -u*dt                                   = ru
                    -c*dt   +A'dy   +  ds   +U'dz           = rd
    -c'dx                   +b'dy           -u'dz   -  dk   = rg
     S*dx                           +X*ds                   = rxs
            +Z*dw                          +W*dz            = rwz
                     k*dt                           +t*dk   = rtk

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
- `rp, rd, ru, rg, rxs, rwz, rtk`: 
"""
function solve_newton_hsd!(
    A, F, b, c, uind::Vector{Int}, uval, θ, θwz,
    p, q, r, ρ,
    x, w, y, s, z, t::RefValue, k::RefValue,
    dx, dw, dy, ds, dz, dt::RefValue, dk::RefValue,
    rp, ru, rd, rg::Float64, rxs, rwz, rtk::Float64
)
    # Solve reduced newton system
    solve_augsys_hsd!(
        A, F, θ, θwz, uind,
        dx, dy, dz,
        rp, 
        rd .- (rxs ./ x),
        ru .- (rwz ./ z)
    )

    # Compute Δτ
    dt.x = (
        rg + (rtk / t.x) + dot(c, dx) - dot(b, dy) + dot(uval, dz)
    ) / ρ
    
    # Compute Δx, Δy, Δz
    dx .+= dt.x .* p
    dy .+= dt.x .* q
    dz .+= dt.x .* r

    # Compute Δs, Δκ, Δw
    # ds  .= (rxs .- s  .* dx)   ./ x
    # dw  .= (rwz .- w  .* dz)   ./ z
    ds  .= rxs
    ds .-= s .* dx
    ds ./= x

    dw  .= rwz
    dw .-= w .* dz
    dw ./= z

    dk.x = (rtk - k.x * dt.x)  / t.x

    return nothing
end

"""
    solve_augsys_hsd!(...)

Solve the augmented system below, and overwrite `dx, dy, dz` with the result.
    -(S/X)*dx  + A'dy  -    U'dz    = rd
         A*dx                       = rp
         U*dx          -(Z/W)*dz    = ru
The augmented system is solved by direct resolution of the normal equations, for
    which a factorization is provided.

# Arguments
- `A`: Matrix
- `F`: Factorization of the matrix `A*Θ*A'`
- `θ`: Diagonal scaling
- `θwz`: Diagonal scaling
- `uind`: Vector of indices of upper-bounded variables
- `dx, dy, dz`: Vectors of unknowns, modified in-place
- `rp, rd, ru`: Right-hand-side vectors
"""
function solve_augsys_hsd!(
    A, F, θ, θwz, uind,
    dx, dy, dz,
    rp, rd, ru
)
    ru_ = ru .* θwz
    rp_ = copy(rd)
    @views rp_[uind] .-= ru_
    rp_ .*= θ
    mul!(dy, A, rp_)
    dy .+= rp

    dy .= (F \ dy)
    mul!(dx, transpose(A), dy)
    dx .-= rd
    @views dx[uind] .+= ru_
    dx .*= θ

    dz .= zero(eltype(dz))
    dz .-= ru
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
function max_step_length(x, dx)
    n = size(x, 1)
    n == size(dx, 1) || throw(DimensionMismatch())
    a = Inf

    @inbounds for i in Base.OneTo(n)
        if dx[i] < 0
            if (-x[i] / dx[i]) < a
                a = (-x[i] / dx[i])
            end
        end
    end
    return a
end

function max_step_length(
    x, w, y, s, z, t::RefValue, k::RefValue,
    dx, dw, dy, ds, dz, dt::RefValue, dk::RefValue
)
    ax = max_step_length(x, dx)
    aw = max_step_length(w, dw)
    as = max_step_length(s, ds)
    az = max_step_length(z, dz)

    at = dt.x < 0.0 ? (-t.x / dt.x) : 1.0
    ak = dk.x < 0.0 ? (-k.x / dk.x) : 1.0
    
    α = min(1.0, ax, aw, as, az, at, ak)

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
    numvar::Int, numvarub::Int, numcon::Int,
    η, γ,
    A, F, b, c, uval, uind::Vector{Int}, θ, θwz,
    p, q, r, ρ,
    x, w, y, s, z, t, k,
    rp, ru, rd, rg, μ,
    dx, dw, dy, ds, dz, dt, dk, α::Float64,
    dxc, dwc, dyc, dsc, dzc, dtc::RefValue, dkc::RefValue,
    beta1::Float64, beta2::Float64, beta3::Float64, beta4::Float64,
)
    # Tentative step length
    α_ = min(1.0, 2.0*α)

    # Tentative cross products
    vx = (x   .+ α_ .* dx)   .* (s   .+ α_ .* ds)
    vw = (w   .+ α_ .* dw)   .* (z   .+ α_ .* dz)
    vt = (t.x  + α_  * dt.x)  * (k.x  + α_  * dk.x)

    # Compute target cross-products
    mu_l = beta4 * μ.x * γ
    mu_u = γ * μ.x / beta4
    @inbounds for i in 1:length(vx)
        if vx[i] < mu_l
            vx[i] = mu_l - vx[i]
        elseif vx[i] > mu_u
            vx[i] = mu_u - vx[i]
        else
            vx[i] = 0.0
        end
    end
    @inbounds for i in 1:length(vw)
        if vw[i] < mu_l
            vw[i] = mu_l - vw[i]
        elseif vx[i] > mu_u
            vw[i] = mu_u - vw[i]
        else
            vw[i] = 0.0
        end
    end
    if vt < mu_l
        vt = mu_l - vt
    elseif vt > mu_u
        vt = mu_u - vt
    else
        vt = 0.0
    end

    # Shift target cross-product to satisfy `v' * e = 0`
    δ = (sum(vx) + sum(vw) + vt) / (numvar + numvarub + 1)
    vx .-= δ
    vw .-= δ
    vt  -= δ

    # Compute corrector
    solve_newton_hsd!(
        A, F, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        x, w, y, s, z, t, k,
        dxc, dwc, dyc, dsc, dzc, dtc, dkc,
        0.0 .* rp, 0.0 .* ru, 0.0 .* rd, 0.0 * rg.x,
        vx,
        vw,
        vt
    )

    # Update corrector
    dxc   .+= dx
    dwc   .+= dw
    dyc   .+= dy
    dsc   .+= ds
    dzc   .+= dz
    dtc.x  += dt.x
    dkc.x  += dk.x

    # Compute corrected step-length
    αc = max_step_length(
        x, w, y, s, z, t, k,
        dxc, dwc, dyc, dsc, dzc, dtc, dkc
    )
    return αc

end