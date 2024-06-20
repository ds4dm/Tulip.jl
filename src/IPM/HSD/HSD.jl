"""
    HSD

Solver for the homogeneous self-dual algorithm.
"""
mutable struct HSD{T, Tv, Tb, Ta, Tk} <: AbstractIPMOptimizer{T}

    # Problem data, in standard form
    dat::IPMData{T, Tv, Tb, Ta}

    # =================
    #   Book-keeping
    # =================
    niter::Int                        # Number of IPM iterations
    solver_status::TerminationStatus  # Optimization status
    primal_status::SolutionStatus     # Primal solution status
    dual_status::SolutionStatus       # Dual   solution status

    primal_objective::T    # Primal objective value: (c'x) / τ
    dual_objective::T      # Dual objective value: (b'y + l' zl - u'zu) / τ

    timer::TimerOutput

    #=====================
        Working memory
    =====================#
    pt::Point{T, Tv}       # Current primal-dual iterate
    res::Residuals{T, Tv}  # Residuals at current iterate
    kkt::Tk
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization
    regG::T   # gap regularization

    function HSD(
        dat::IPMData{T, Tv, Tb, Ta}, kkt_options::KKTOptions{T}
    ) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}}

        m, n = dat.nrow, dat.ncol
        p = sum(dat.lflag) + sum(dat.uflag)

        # Allocate some memory
        pt  = Point{T, Tv}(m, n, p, hflag=true)
        res = Residuals(
            tzeros(Tv, m), tzeros(Tv, n), tzeros(Tv, n),
            tzeros(Tv, n), zero(T),
            zero(T), zero(T), zero(T), zero(T), zero(T)
        )

        # Initial regularizations
        regP = tones(Tv, n)
        regD = tones(Tv, m)
        regG = one(T)

        kkt = KKT.setup(dat.A, kkt_options.System, kkt_options.Backend)
        Tk = typeof(kkt)

        return new{T, Tv, Tb, Ta, Tk}(dat,
            0, Trm_Unknown, Sln_Unknown, Sln_Unknown,
            T(Inf), T(-Inf),
            TimerOutput(),
            pt, res, kkt, regP, regD, regG
        )
    end

end

include("dot_for_mutable.jl")

include("step.jl")


"""
    compute_residuals!(::HSD, res, pt, A, b, c, uind, uval)

In-place computation of primal-dual residuals at point `pt`.
"""
# TODO: check whether having just hsd as argument makes things slower
# TODO: Update solution status
function compute_residuals!(hsd::HSD{T}
) where{T}

    pt, res = hsd.pt, hsd.res
    dat = hsd.dat

    # Primal residual
    # rp = t*b - A*x
    axpby!(pt.τ, dat.b, zero(T), res.rp)
    mul!(res.rp, dat.A, pt.x, -one(T), one(T))

    # Lower-bound residual
    # rl_j = τ*l_j - (x_j - xl_j)  if l_j ∈ R
    #      = 0                     if l_j = -∞
    @. res.rl = (- pt.x + pt.xl + pt.τ * dat.l) * dat.lflag

    # Upper-bound residual
    # ru_j = τ*u_j - (x_j + xu_j)  if u_j ∈ R
    #      = 0                     if u_j = +∞
    @. res.ru = (- pt.x - pt.xu + pt.τ * dat.u) * dat.uflag

    # Dual residual
    # rd = t*c - A'y - zl + zu
    axpby!(pt.τ, dat.c, zero(T), res.rd)
    mul!(res.rd, transpose(dat.A), pt.y, -one(T), one(T))
    @. res.rd += pt.zu .* dat.uflag - pt.zl .* dat.lflag

    dot_buf = buffer_for_dot_weighted_sum(T)

    # Gap residual
    # rg = c'x - (b'y + l'zl - u'zu) + k
    res.rg = pt.κ + buffered_dot_weighted_sum!!(
        dot_buf,
        (
            (dat.c, pt.x),
            (dat.b, pt.y),
            (dat.l[dat.lflag], pt.zl[dat.lflag]),
            (dat.u[dat.uflag], pt.zu[dat.uflag]),
        ),
        (
            1, -1, -1, 1,
        ),
    )

    # Residuals norm
    res.rp_nrm = norm(res.rp, Inf)
    res.rl_nrm = norm(res.rl, Inf)
    res.ru_nrm = norm(res.ru, Inf)
    res.rd_nrm = norm(res.rd, Inf)
    res.rg_nrm = norm(res.rg, Inf)

    # Compute primal and dual bounds
    hsd.primal_objective = buffered_dot_product!!(dot_buf.dot, dat.c, pt.x) / pt.τ + dat.c0
    hsd.dual_objective = buffered_dot_weighted_sum!!(
        dot_buf,
        (
            (dat.b, pt.y),
            (dat.l[dat.lflag], pt.zl[dat.lflag]),
            (dat.u[dat.uflag], pt.zu[dat.uflag]),
        ),
        (
            1, 1, -1,
        ),
    ) / pt.τ + dat.c0

    return nothing
end


"""
    update_solver_status!()

Update status and return true if solver should stop.
"""
function update_solver_status!(hsd::HSD{T}, ϵp::T, ϵd::T, ϵg::T, ϵi::T) where{T}
    hsd.solver_status = Trm_Unknown

    pt, res = hsd.pt, hsd.res
    dat = hsd.dat

    ρp = max(
        res.rp_nrm / (pt.τ * (one(T) + norm(dat.b, Inf))),
        res.rl_nrm / (pt.τ * (one(T) + norm(dat.l .* dat.lflag, Inf))),
        res.ru_nrm / (pt.τ * (one(T) + norm(dat.u .* dat.uflag, Inf)))
    )
    ρd = res.rd_nrm / (pt.τ * (one(T) + norm(dat.c, Inf)))
    ρg = abs(hsd.primal_objective - hsd.dual_objective) / (one(T) + abs(hsd.dual_objective))

    # Check for feasibility
    if ρp <= ϵp
        hsd.primal_status = Sln_FeasiblePoint
    else
        hsd.primal_status = Sln_Unknown
    end

    if ρd <= ϵd
        hsd.dual_status = Sln_FeasiblePoint
    else
        hsd.dual_status = Sln_Unknown
    end

    # Check for optimal solution
    if ρp <= ϵp && ρd <= ϵd && ρg <= ϵg
        hsd.primal_status = Sln_Optimal
        hsd.dual_status   = Sln_Optimal
        hsd.solver_status = Trm_Optimal
        return nothing
    end

    dot_buf = buffer_for_dot_weighted_sum(T)

    # Check for infeasibility certificates
    if max(
        norm(dat.A * pt.x, Inf),
        norm((pt.x .- pt.xl) .* dat.lflag, Inf),
        norm((pt.x .+ pt.xu) .* dat.uflag, Inf)
    ) * (norm(dat.c, Inf) / max(1, norm(dat.b, Inf))) <
        -ϵi * buffered_dot_product!!(dot_buf.dot, dat.c, pt.x)
        # Dual infeasible, i.e., primal unbounded
        hsd.primal_status = Sln_InfeasibilityCertificate
        hsd.solver_status = Trm_DualInfeasible
        return nothing
    end

    δ = dat.A' * pt.y .+ (pt.zl .* dat.lflag) .- (pt.zu .* dat.uflag)
    if norm(δ, Inf) * max(
        norm(dat.l .* dat.lflag, Inf),
        norm(dat.u .* dat.uflag, Inf),
        norm(dat.b, Inf)
    ) / (max(one(T), norm(dat.c, Inf)))  < buffered_dot_weighted_sum!!(
        dot_buf,
        (
            (dat.b, pt.y),
            (dat.l[dat.lflag], pt.zl[dat.lflag]),
            (dat.u[dat.uflag], pt.zu[dat.uflag]),
        ),
        (
            1, 1, -1,
        ),
    ) * ϵi

        # Primal infeasible
        hsd.dual_status = Sln_InfeasibilityCertificate
        hsd.solver_status = Trm_PrimalInfeasible
        return nothing
    end

    return nothing
end


"""
    optimize!

"""
function ipm_optimize!(hsd::HSD{T}, params::IPMOptions{T}) where{T}
    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function
    dat = hsd.dat

    # Initialization
    TimerOutputs.reset_timer!(hsd.timer)
    tstart = time()
    hsd.niter = 0

    # Print information about the problem
    if params.OutputLevel > 0
        @printf "\nOptimizer info (HSD)\n"
        @printf "Constraints  : %d\n" dat.nrow
        @printf "Variables    : %d\n" dat.ncol
        bmin, bmax = extrema(dat.b)
        @printf "RHS          : [%+.2e, %+.2e]\n" bmin bmax
        lmin, lmax = extrema(dat.l .* dat.lflag)
        @printf "Lower bounds : [%+.2e, %+.2e]\n" lmin lmax
        lmin, lmax = extrema(dat.u .* dat.uflag)
        @printf "Upper bounds : [%+.2e, %+.2e]\n" lmin lmax


        @printf "\nLinear solver options\n"
        @printf "  %-12s : %s\n" "Arithmetic" KKT.arithmetic(hsd.kkt)
        @printf "  %-12s : %s\n" "Backend" KKT.backend(hsd.kkt)
        @printf "  %-12s : %s\n" "System" KKT.linear_system(hsd.kkt)
    end

    # IPM LOG
    if params.OutputLevel > 0
        @printf "\n%4s  %14s  %14s  %8s %8s %8s  %7s  %4s\n" "Itn" "PObj" "DObj" "PFeas" "DFeas" "GFeas" "Mu" "Time"
    end

    # Set starting point
    hsd.pt.x   .= zero(T)
    hsd.pt.xl  .= one(T) .* dat.lflag
    hsd.pt.xu  .= one(T) .* dat.uflag

    hsd.pt.y   .= zero(T)
    hsd.pt.zl  .= one(T) .* dat.lflag
    hsd.pt.zu  .= one(T) .* dat.uflag

    hsd.pt.τ   = one(T)
    hsd.pt.κ   = one(T)

    update_mu!(hsd.pt)

    # Main loop
    # Iteration 0 corresponds to the starting point.
    # Therefore, there is no numerical factorization before the first log is printed.
    # If the maximum number of iterations is set to 0, the only computation that occurs
    # is computing the residuals at the initial point.
    @timeit hsd.timer "Main loop" while(true)

        # I.A - Compute residuals at current iterate
        @timeit hsd.timer "Residuals" compute_residuals!(hsd)

        update_mu!(hsd.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if params.OutputLevel > 0
            # Display log
            @printf "%4d" hsd.niter

            # Objectives
            ϵ = dat.objsense ? one(T) : -one(T)
            @printf "  %+14.7e" ϵ * hsd.primal_objective
            @printf "  %+14.7e" ϵ * hsd.dual_objective

            # Residuals
            @printf "  %8.2e" max(hsd.res.rp_nrm, hsd.res.ru_nrm)
            @printf " %8.2e" hsd.res.rd_nrm
            @printf " %8.2e" hsd.res.rg_nrm

            # Mu
            @printf "  %7.1e" hsd.pt.μ

            # Time
            @printf "  %.2f" ttot

            print("\n")
        end

        # TODO: check convergence status
        # TODO: first call an `compute_convergence status`,
        #   followed by a check on the solver status to determine whether to stop
        # In particular, user limits should be checked last (if an optimal solution is found,
        # we want to report optimal, not user limits)
        @timeit hsd.timer "update status" update_solver_status!(hsd,
            params.TolerancePFeas,
            params.ToleranceDFeas,
            params.ToleranceRGap,
            params.ToleranceIFeas
        )

        if (
            hsd.solver_status == Trm_Optimal
            || hsd.solver_status == Trm_PrimalInfeasible
            || hsd.solver_status == Trm_DualInfeasible
        )
            break
        elseif hsd.niter >= params.IterationsLimit
            hsd.solver_status = Trm_IterationLimit
            break
        elseif ttot >= params.TimeLimit
            hsd.solver_status = Trm_TimeLimit
            break
        end


        # TODO: step
        # For now, include the factorization in the step function
        # Q: should we use more arguments here?
        try
            @timeit hsd.timer "Step" compute_step!(hsd, params)
        catch err

            if isa(err, PosDefException) || isa(err, SingularException)
                # Numerical trouble while computing the factorization
                hsd.solver_status = Trm_NumericalProblem

            elseif isa(err, OutOfMemoryError)
                # Out of memory
                hsd.solver_status = Trm_MemoryLimit

            elseif isa(err, InterruptException)
                hsd.solver_status = Trm_Unknown
            else
                # Unknown error: rethrow
                rethrow(err)
            end

            break
        end

        hsd.niter += 1

    end

    # TODO: print message based on termination status
    params.OutputLevel > 0 && println("Solver exited with status $((hsd.solver_status))")

    return nothing

end
