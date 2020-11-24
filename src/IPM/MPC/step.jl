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
    mpc.regP ./= 10
    mpc.regD ./= 10
    clamp!(mpc.regP, sqrt(eps(T)), one(T))
    clamp!(mpc.regD, sqrt(eps(T)), one(T))

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

    # II. Compute search direction
    Δ = mpc.Δ
    
    # Affine-scaling direction (predictor)
    @timeit mpc.timer "Newton" solve_newton_system!(Δ, mpc,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd,
        -(pt.xl .* pt.zl) .* dat.lflag,
        -(pt.xu .* pt.zu) .* dat.uflag,
    )

    # Step length for affine-scaling direction
    αp, αd = max_step_length_pd(pt, Δ)
    μₐ = (
        dot((pt.xl + αp .* Δ.xl) .* dat.lflag, pt.zl + αd .* Δ.zl)
        + dot((pt.xu + αp .* Δ.xu) .* dat.uflag, pt.zu + αd .* Δ.zu)
    ) / pt.p
    σ = clamp((μₐ / pt.μ)^3, sqrt(eps(T)), one(T) - sqrt(eps(T)))
    
    # Mehrotra corrector
    @timeit mpc.timer "Newton" solve_newton_system!(Δ, mpc,
        # Right-hand side of Newton system
        res.rp, res.rl, res.ru, res.rd,
        (-(pt.xl .* pt.zl) .+ σ * pt.μ .- Δ.xl .* Δ.zl) .* dat.lflag,
        (-(pt.xu .* pt.zu) .+ σ * pt.μ .- Δ.xu .* Δ.zu) .* dat.uflag,
    )
    αp, αd = max_step_length_pd(pt, Δ)

    # Update current iterate
    αp *= params.StepDampFactor
    αd *= params.StepDampFactor
    pt.x  .+= αp .* Δ.x
    pt.xl .+= αp .* Δ.xl
    pt.xu .+= αp .* Δ.xu
    pt.y  .+= αd .* Δ.y
    pt.zl .+= αd .* Δ.zl
    pt.zu .+= αd .* Δ.zu
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

"""
    max_step_length_pd(pt, δ)

Compute maximum primal-dual step length.
"""
function max_step_length_pd(pt::Point{T, Tv}, δ::Point{T, Tv}) where{T, Tv<:AbstractVector{T}}
    axl = max_step_length(pt.xl, δ.xl)
    axu = max_step_length(pt.xu, δ.xu)
    azl = max_step_length(pt.zl, δ.zl)
    azu = max_step_length(pt.zu, δ.zu)

    αp = min(one(T), axl, axu)
    αd = min(one(T), azl, azu)

    return αp, αd
end