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
    Δ  = mpc.Δ
    Δc = mpc.Δc

    # Affine-scaling direction and associated step size
    @timeit mpc.timer "Predictor" compute_predictor!(mpc::MPC)
    mpc.αp, mpc.αd = max_step_length_pd(mpc.pt, mpc.Δ)

    # TODO: if step size is large enough, skip corrector

    # Corrector
    @timeit mpc.timer "Corrector" compute_corrector!(mpc::MPC)
    mpc.αp, mpc.αd = max_step_length_pd(mpc.pt, mpc.Δc)

    # Update current iterate
    mpc.αp *= params.StepDampFactor
    mpc.αd *= params.StepDampFactor
    pt.x  .+= mpc.αp .* Δc.x
    pt.xl .+= mpc.αp .* Δc.xl
    pt.xu .+= mpc.αp .* Δc.xu
    pt.y  .+= mpc.αd .* Δc.y
    pt.zl .+= mpc.αd .* Δc.zl
    pt.zu .+= mpc.αd .* Δc.zu
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
    # @printf "|rp|   = %16.8e\n" norm(dat.A * Δ.x - ξp, Inf)
    # @printf "|rl|   = %16.8e\n" norm((Δ.x - Δ.xl) .* dat.lflag - ξl, Inf)
    # @printf "|ru|   = %16.8e\n" norm((Δ.x + Δ.xu) .* dat.uflag - ξu, Inf)
    # @printf "|rd|   = %16.8e\n" norm(dat.A'Δ.y + Δ.zl - Δ.zu - ξd, Inf)
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


"""
    compute_predictor!(mpc::MPC) -> Nothing
"""
function compute_predictor!(mpc::MPC)

    # Newton RHS
    copyto!(mpc.ξp, mpc.res.rp)
    copyto!(mpc.ξl, mpc.res.rl)
    copyto!(mpc.ξu, mpc.res.ru)
    copyto!(mpc.ξd, mpc.res.rd)
    @. mpc.ξxzl = -(mpc.pt.xl .* mpc.pt.zl) .* mpc.dat.lflag
    @. mpc.ξxzu = -(mpc.pt.xu .* mpc.pt.zu) .* mpc.dat.uflag

    # Compute affine-scaling direction
    @timeit mpc.timer "Newton" solve_newton_system!(mpc.Δ, mpc,
        mpc.ξp, mpc.ξl, mpc.ξu, mpc.ξd, mpc.ξxzl, mpc.ξxzu
    )

    # TODO: check Newton system residuals, perform iterative refinement if needed
    return nothing
end

"""
    compute_corrector!(mpc::MPC) -> Nothing
"""
function compute_corrector!(mpc::MPC{T, Tv}) where{T, Tv<:AbstractVector{T}}
    dat = mpc.dat
    pt = mpc.pt
    Δ = mpc.Δ
    Δc = mpc.Δc

    # Step length for affine-scaling direction
    αp_aff, αd_aff = mpc.αp, mpc.αd
    μₐ = (
        dot((pt.xl + αp_aff .* Δ.xl) .* dat.lflag, pt.zl + αd_aff .* Δ.zl)
        + dot((pt.xu + αp_aff .* Δ.xu) .* dat.uflag, pt.zu + αd_aff .* Δ.zu)
    ) / pt.p
    σ = clamp((μₐ / pt.μ)^3, sqrt(eps(T)), one(T) - sqrt(eps(T)))

    # Newton RHS
    # rmul!(mpc.ξp, zero(T))
    # rmul!(mpc.ξl, zero(T))
    # rmul!(mpc.ξu, zero(T))
    # rmul!(mpc.ξd, zero(T))
    @. mpc.ξxzl = (σ * pt.μ .- Δ.xl .* Δ.zl .- pt.xl .* pt.zl) .* dat.lflag
    @. mpc.ξxzu = (σ * pt.μ .- Δ.xu .* Δ.zu .- pt.xu .* pt.zu) .* dat.uflag

    # Compute corrector
    @timeit mpc.timer "Newton" solve_newton_system!(mpc.Δc, mpc,
        mpc.ξp, mpc.ξl, mpc.ξu, mpc.ξd, mpc.ξxzl, mpc.ξxzu
    )

    # TODO: check Newton system residuals, perform iterative refinement if needed
    return nothing
end
