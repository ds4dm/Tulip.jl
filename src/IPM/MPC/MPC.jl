"""
    MPC

Implements Mehrotra's Predictor-Corrector interior-point algorithm.
"""
mutable struct MPC{T, Tv, Tb, Ta, Tk} <: AbstractIPMOptimizer{T}

    # Problem data, in standard form
    dat::IPMData{T, Tv, Tb, Ta}

    # =================
    #   Book-keeping
    # =================
    niter::Int                        # Number of IPM iterations
    solver_status::TerminationStatus  # Optimization status
    primal_status::SolutionStatus     # Primal solution status
    dual_status::SolutionStatus       # Dual   solution status


    primal_objective::T  # Primal bound: c'x
    dual_objective::T    # Dual bound: b'y + l' zl - u'zu

    timer::TimerOutput

    #=====================
        Working memory
    =====================#
    pt::Point{T, Tv}       # Current primal-dual iterate
    res::Residuals{T, Tv}  # Residuals at current iterate

    Δ::Point{T, Tv}   # Predictor direction
    Δc::Point{T, Tv}  # Corrector

    # Step sizes
    αp::T
    αd::T

    # Newton system RHS

    kkt::Tk
    regP::Tv  # Primal regularization
    regD::Tv  # Dual regularization

    function MPC(
        dat::IPMData{T, Tv, Tb, Ta}, kkt_options::KKTOptions{T}
    ) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}}
        
        m, n = dat.nrow, dat.ncol
        p = sum(dat.lflag) + sum(dat.uflag)

        # Allocate some memory
        pt  = Point{T, Tv}(m, n, p, hflag=false)
        res = Residuals(
            tzeros(Tv, m), tzeros(Tv, n), tzeros(Tv, n),
            tzeros(Tv, n), zero(T),
            zero(T), zero(T), zero(T), zero(T), zero(T)
        )
        Δ  = Point{T, Tv}(m, n, p, hflag=false)
        Δc = Point{T, Tv}(m, n, p, hflag=false)

        # Initial regularizations
        regP = tones(Tv, n)
        regD = tones(Tv, m)

        kkt = KKT.setup(kkt_options.Factory.T, dat.A; kkt_options.Factory.options...)
        Tk = typeof(kkt)

        return new{T, Tv, Tb, Ta, Tk}(dat,
            0, Trm_Unknown, Sln_Unknown, Sln_Unknown,
            T(Inf), T(-Inf),
            TimerOutput(),
            pt, res, Δ, Δc, zero(T), zero(T),
            kkt, regP, regD
        )
    end

end

include("step.jl")

"""
    compute_residuals!(::MPC)

In-place computation of primal-dual residuals at point `pt`.
"""
function compute_residuals!(mpc::MPC{T}) where{T}

    pt, res = mpc.pt, mpc.res
    dat = mpc.dat

    # Primal residual
    # rp = b - A*x
    res.rp .= dat.b
    mul!(res.rp, dat.A, pt.x, -one(T), one(T))

    # Lower-bound residual
    # rl_j = l_j - (x_j - xl_j)  if l_j ∈ R
    #      = 0                   if l_j = -∞
    @. res.rl = (dat.l + pt.xl) - pt.x
    res.rl .*= dat.lflag

    # Upper-bound residual
    # ru_j = u_j - (x_j + xu_j)  if u_j ∈ R
    #      = 0                   if u_j = +∞
    @. res.ru = dat.u - (pt.x + pt.xu)
    res.ru .*= dat.uflag

    # Dual residual
    # rd = c - (A'y + zl - zu)
    res.rd .= dat.c
    mul!(res.rd, transpose(dat.A), pt.y, -one(T), one(T))
    @. res.rd += pt.zu .* dat.uflag - pt.zl .* dat.lflag

    # Residuals norm
    res.rp_nrm = norm(res.rp, Inf)
    res.rl_nrm = norm(res.ru, Inf)
    res.ru_nrm = norm(res.ru, Inf)
    res.rd_nrm = norm(res.rd, Inf)

    # Compute primal and dual bounds
    mpc.primal_objective = dot(dat.c, pt.x) + dat.c0
    mpc.dual_objective   = (
        dot(dat.b, pt.y)
        + dot(dat.l .* dat.lflag, pt.zl)
        - dot(dat.u .* dat.uflag, pt.zu)
    ) + dat.c0

    return nothing
end


"""
    update_solver_status!()

Update status and return true if solver should stop.
"""
function update_solver_status!(mpc::MPC{T}, ϵp::T, ϵd::T, ϵg::T, ϵi::T) where{T}
    mpc.solver_status = Trm_Unknown

    pt, res = mpc.pt, mpc.res
    dat = mpc.dat

    ρp = max(
        res.rp_nrm / (one(T) + norm(dat.b, Inf)),
        res.rl_nrm / (one(T) + norm(dat.l .* dat.lflag, Inf)),
        res.ru_nrm / (one(T) + norm(dat.u .* dat.uflag, Inf))
    )
    ρd = res.rd_nrm / (one(T) + norm(dat.c, Inf))
    ρg = abs(mpc.primal_objective - mpc.dual_objective) / (one(T) + abs(mpc.primal_objective))

    # Check for feasibility
    if ρp <= ϵp
        mpc.primal_status = Sln_FeasiblePoint
    else
        mpc.primal_status = Sln_Unknown
    end

    if ρd <= ϵd
        mpc.dual_status = Sln_FeasiblePoint
    else
        mpc.dual_status = Sln_Unknown
    end
    
    # Check for optimal solution
    if ρp <= ϵp && ρd <= ϵd && ρg <= ϵg
        mpc.primal_status = Sln_Optimal
        mpc.dual_status   = Sln_Optimal
        mpc.solver_status = Trm_Optimal
        return nothing
    end
    
    # TODO: Primal/Dual infeasibility detection
    # Check for infeasibility certificates
    if max(
        norm(dat.A * pt.x, Inf),
        norm((pt.x - pt.xl) .* dat.lflag, Inf),
        norm((pt.x + pt.xu) .* dat.uflag, Inf)
    ) * (norm(dat.c, Inf) / max(1, norm(dat.b, Inf))) < - ϵi * dot(dat.c, pt.x)
        # Dual infeasible, i.e., primal unbounded
        mpc.primal_status = Sln_InfeasibilityCertificate
        mpc.solver_status = Trm_DualInfeasible
        return nothing
    end

    δ = dat.A' * pt.y + (pt.zl .* dat.lflag) - (pt.zu .* dat.uflag)
    if norm(δ, Inf) * max(
        norm(dat.l .* dat.lflag, Inf),
        norm(dat.u .* dat.uflag, Inf),
        norm(dat.b, Inf)
    ) / (max(one(T), norm(dat.c, Inf)))  < (dot(dat.b, pt.y) + dot(dat.l .* dat.lflag, pt.zl)- dot(dat.u .* dat.uflag, pt.zu)) * ϵi
        # Primal infeasible
        mpc.dual_status = Sln_InfeasibilityCertificate
        mpc.solver_status = Trm_PrimalInfeasible
        return nothing
    end

    return nothing
end


"""
    optimize!

"""
function ipm_optimize!(mpc::MPC{T}, params::IPMOptions{T}) where{T}
    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function
    dat = mpc.dat

    # Initialization
    TimerOutputs.reset_timer!(mpc.timer)
    tstart = time()
    mpc.niter = 0

    # Print information about the problem
    if params.OutputLevel > 0
        @printf "\nOptimizer info (MPC)\n"
        @printf "Constraints  : %d\n" dat.nrow
        @printf "Variables    : %d\n" dat.ncol
        bmin, bmax = extrema(dat.b)
        @printf "RHS          : [%+.2e, %+.2e]\n" bmin bmax
        lmin, lmax = extrema(dat.l .* dat.lflag)
        @printf "Lower bounds : [%+.2e, %+.2e]\n" lmin lmax
        lmin, lmax = extrema(dat.u .* dat.uflag)
        @printf "Upper bounds : [%+.2e, %+.2e]\n" lmin lmax


        @printf "\nLinear solver options\n"
        @printf "  %-12s : %s\n" "Arithmetic" KKT.arithmetic(mpc.kkt)
        @printf "  %-12s : %s\n" "Backend" KKT.backend(mpc.kkt)
        @printf "  %-12s : %s\n" "System" KKT.linear_system(mpc.kkt)
    end

    # IPM LOG
    if params.OutputLevel > 0
        @printf "\n%4s  %14s  %14s  %8s %8s %8s  %7s  %4s\n" "Itn" "PObj" "DObj" "PFeas" "DFeas" "GFeas" "Mu" "Time"
    end

    # Set starting point
    @timeit mpc.timer "Initial point" compute_starting_point(mpc)

    # Main loop
    # Iteration 0 corresponds to the starting point.
    # Therefore, there is no numerical factorization before the first log is printed.
    # If the maximum number of iterations is set to 0, the only computation that occurs
    # is computing the residuals at the initial point.
    @timeit mpc.timer "Main loop" while(true)

        # I.A - Compute residuals at current iterate
        @timeit mpc.timer "Residuals" compute_residuals!(mpc)

        update_mu!(mpc.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if params.OutputLevel > 0
            # Display log
            @printf "%4d" mpc.niter
            
            # Objectives
            ϵ = dat.objsense ? one(T) : -one(T)
            @printf "  %+14.7e" ϵ * mpc.primal_objective
            @printf "  %+14.7e" ϵ * mpc.dual_objective
            
            # Residuals
            @printf "  %8.2e" max(mpc.res.rp_nrm, mpc.res.rl_nrm, mpc.res.ru_nrm)
            @printf " %8.2e" mpc.res.rd_nrm
            @printf " %8s" "--"

            # Mu
            @printf "  %7.1e" mpc.pt.μ

            # Time
            @printf "  %.2f" ttot

            print("\n")
        end

        # TODO: check convergence status
        # TODO: first call an `compute_convergence status`,
        #   followed by a check on the solver status to determine whether to stop
        # In particular, user limits should be checked last (if an optimal solution is found,
        # we want to report optimal, not user limits)
        @timeit mpc.timer "update status" update_solver_status!(mpc,
            params.BarrierTolerancePFeas,
            params.BarrierToleranceDFeas,
            params.BarrierToleranceRGap,
            params.BarrierToleranceIFeas
        )

        if (
            mpc.solver_status == Trm_Optimal
            || mpc.solver_status == Trm_PrimalInfeasible
            || mpc.solver_status == Trm_DualInfeasible
        )
            break
        elseif mpc.niter >= params.BarrierIterationsLimit 
            mpc.solver_status = Trm_IterationLimit
            break
        elseif ttot >= params.TimeLimit
            mpc.solver_status = Trm_TimeLimit
            break
        end
        

        # TODO: step
        # For now, include the factorization in the step function
        # Q: should we use more arguments here?
        try
            @timeit mpc.timer "Step" compute_step!(mpc, params)
        catch err

            if isa(err, PosDefException) || isa(err, SingularException)
                # Numerical trouble while computing the factorization
                mpc.solver_status = Trm_NumericalProblem
    
            elseif isa(err, OutOfMemoryError)
                # Out of memory
                mpc.solver_status = Trm_MemoryLimit

            elseif isa(err, InterruptException)
                mpc.solver_status = Trm_Unknown
            else
                # Unknown error: rethrow
                rethrow(err)
            end

            break
        end

        mpc.niter += 1

    end
    
    # TODO: print message based on termination status
    params.OutputLevel > 0 && println("Solver exited with status $((mpc.solver_status))")

    return nothing

end

function compute_starting_point(mpc::MPC{T}) where{T}

    pt = mpc.pt
    dat = mpc.dat
    m, n, p = pt.m, pt.n, pt.p

    KKT.update!(mpc.kkt, zeros(T, n), ones(T, n), T(1e-6) .* ones(T, m))

    # Get initial iterate
    KKT.solve!(zeros(T, n), pt.y, mpc.kkt, false .* mpc.dat.b, mpc.dat.c)  # For y
    KKT.solve!(pt.x, zeros(T, m), mpc.kkt, mpc.dat.b, false .* mpc.dat.c)  # For x
    
    # I. Recover positive primal-dual coordinates
    δx = one(T) + max(
        zero(T),
        (-3 // 2) * minimum((pt.x - dat.l) .* dat.lflag),
        (-3 // 2) * minimum((dat.u - pt.x) .* dat.uflag)
    )
    pt.xl  .= ((pt.x - dat.l) .+ δx) .* dat.lflag
    pt.xu  .= ((dat.u - pt.x) .+ δx) .* dat.uflag

    z = dat.c - dat.A' * pt.y
    #=
        We set zl, zu such that `z = zl - zu`

         lⱼ |  uⱼ |    zˡⱼ |     zᵘⱼ |
        ----+-----+--------+---------+
        yes | yes | ¹/₂ zⱼ | ⁻¹/₂ zⱼ |
        yes |  no |     zⱼ |      0  |
         no | yes |     0  |     -zⱼ |
         no |  no |     0  |      0  |
        ----+-----+--------+---------+
    =#
    pt.zl .= ( z ./ (dat.lflag + dat.uflag)) .* dat.lflag
    pt.zu .= (-z ./ (dat.lflag + dat.uflag)) .* dat.uflag
    
    δz = one(T) + max(zero(T), (-3 // 2) * minimum(pt.zl), (-3 // 2) * minimum(pt.zu))
    pt.zl[dat.lflag] .+= δz
    pt.zu[dat.uflag] .+= δz
    
    mpc.pt.τ   = one(T)
    mpc.pt.κ   = zero(T)

    # II. Balance complementarity products
    μ = dot(pt.xl, pt.zl) + dot(pt.xu, pt.zu)
    dx = μ / ( 2 * (sum(pt.zl) + sum(pt.zu)))
    dz = μ / ( 2 * (sum(pt.xl) + sum(pt.xu)))

    pt.xl[dat.lflag] .+= dx
    pt.xu[dat.uflag] .+= dx
    pt.zl[dat.lflag] .+= dz
    pt.zu[dat.uflag] .+= dz

    # Update centrality parameter
    update_mu!(mpc.pt)

    return nothing
end