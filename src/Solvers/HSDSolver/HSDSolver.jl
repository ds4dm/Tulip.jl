"""
    HSDSolver

Solver for the homogeneous self-dual algorithm.
"""
mutable struct HSDSolver{T, Tv, Tk} <: AbstractIPMSolver{T}

    # =================
    #   Book-keeping
    # =================
    niter::Int                        # Number of IPM iterations
    solver_status::TerminationStatus  # Optimization status
    primal_status::SolutionStatus     # Primal solution status
    dual_status::SolutionStatus       # Dual   solution status

    primal_bound_unscaled::T  # Unscaled primal bound c'x
    primal_bound_scaled::T    # Scaled primal bound (c'x) / τ
    dual_bound_unscaled::T    # Unscaled dual bound b'y + l' zl - u'zu
    dual_bound_scaled::T      # Scaled dual bound (b'y + l' zl - u'zu) / τ

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

    function HSDSolver(
        dat::IPMData{T, Tv, Tb, Ta}, params::Parameters{T}
    ) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}}
        
        m, n = dat.nrow, dat.ncol
        p = sum(dat.lflag) + sum(dat.uflag)

        # Allocate some memory
        pt  = Point{T, Tv}(m, n, p)
        res = Residuals(
            tzeros(Tv, m), tzeros(Tv, n), tzeros(Tv, n),
            tzeros(Tv, n), zero(T),
            zero(T), zero(T), zero(T), zero(T), zero(T)
        )

        # Initial regularizations
        regP = tones(Tv, n)
        regD = tones(Tv, m)
        regG = one(T)

        kkt = KKT.setup(params.KKTOptions.Ts, dat.A; params.KKTOptions.options...)
        Tk = typeof(kkt)

        return new{T, Tv, Tk}(
            0, Trm_Unknown, Sln_Unknown, Sln_Unknown,
            T(Inf), T(Inf), T(-Inf), T(-Inf),
            TimerOutput(),
            pt, res, kkt, regP, regD, regG
        )
    end

    # TODO: remove and convert ProblemData --> IPMData
    function HSDSolver{Tv}(params::Parameters{Tv}, pb::ProblemData{Tv}) where{Tv}
        nvar = pb.nvar
        ncon = pb.ncon

        # =============================================
        #   Convert problem to standard form
        # =============================================

        # Count memory needed
        nz = 0  # Non-zeros (excluding slacks)
        nfree = 0  # Number of free variables
        nvarupb = 0  # Upper-bounded variables (excluding slacks)
        nslack = 0  # Slack variables
        nslackupb = 0  # Upper-bounded slacks

        # Count slacks
        for (i, (l, u)) in enumerate(zip(pb.lcon, pb.ucon))
            if l == u
                # a'x == b
            elseif l == -Inf && isfinite(u)
                # `a'x ⩽ u` becomes `a'x + s == u, s ⩾ 0`
                nslack += 1
            elseif isfinite(l) && u == Inf
                # `a'x ⩾ l` becomes `a'x - s == l, s ⩾ 0`
                nslack += 1
            elseif isfinite(l) && isfinite(u)
                # `l ⩽ a'x ⩽ u` becomes `a'x - s == l, 0 ⩽ s ⩽ u-l`
                nslack += 1
                nslackupb += 1
            else
                error("Invalid bounds ($l, $u) for row $i")
            end
        end

        # Count non-zeros coeffs, free & upper-bounded variables,
        # and flag which columns need to be flipped
        # TODO: handle free variables explicitly
        nz = 0
        for (j, (l, u)) in enumerate(zip(pb.lvar, pb.uvar))
            if isfinite(l) && isfinite(u)
                # l ⩽ x ⩽ u
                nvarupb += 1
            elseif l == -Inf && isfinite(u)
                # x ⩽ u
            elseif isfinite(l) && u == Inf
                # l ⩽ x
            elseif l == -Inf && u == Inf
                # Free variable
                # x is replaced by `x⁺ - x⁻`
                nfree += 1
                nz += length(pb.acols[j].nzind)
            else
                error("Invalid bounds ($l, $u) for variable $j")
            end
            nz += length(pb.acols[j].nzind)
        end

        # Allocate memory
        c = Vector{Tv}(undef, nvar + nfree + nslack)
        c0 = pb.obj0
        b = Vector{Tv}(undef, ncon)
        aI = Vector{Int}(undef, nz + nslack)
        aJ = Vector{Int}(undef, nz + nslack)
        aV = Vector{Tv}(undef, nz + nslack)
        uind = Vector{Int}(undef, nvarupb + nslackupb)
        uval = Vector{Tv}(undef, nvarupb + nslackupb)

        # Populate right-hand side, slack coefficients and bounds
        # This needs to happen before we populate the matrix,
        # otherwise the right-hand side values are meaningless
        jslack = nvar + nfree  # index of slack variable
        jslackupb = nvarupb
        # @info "Non-zeros before adding slacks" nz
        for (i, (l, u)) in enumerate(zip(pb.lcon, pb.ucon))
            if l == u
                # `a'x == b`, nothing to do
                b[i] = l
            elseif l == -Inf && isfinite(u)
                # @info "Row: a'x <= $u"
                # `a'x ⩽ u` becomes `a'x + s == u, s ⩾ 0`
                nz += 1
                jslack += 1
                # Right-hand side
                b[i] = u
                # Slack coefficient
                aI[nz] = i
                aJ[nz] = jslack
                aV[nz] = one(Tv)
                # Slack objective
                c[jslack] = zero(Tv)
            elseif isfinite(l) && u == Inf
                # `a'x ⩾ l` becomes `a'x - s == l, s ⩾ 0`
                nz += 1
                jslack += 1
                # Right-hand side
                b[i] = l
                # Slack coefficient
                aI[nz] = i
                aJ[nz] = jslack
                aV[nz] = -one(Tv)
                # Slack objective
                c[jslack] = zero(Tv)
            elseif isfinite(l) && isfinite(u)
                # `l ⩽ a'x ⩽ u` becomes `a'x - s == l, 0 ⩽ s ⩽ u-l`
                nz += 1
                jslack += 1
                jslackupb += 1
                # Right-hand side
                b[i] = l
                # Slack coefficient
                aI[nz] = i
                aJ[nz] = jslack
                aV[nz] = -one(Tv)
                # Slack objective
                c[jslack] = zero(Tv)
                # TODO: slack upper bound
                uind[jslackupb] = jslack
                uval[jslackupb] = u - l
            else
                error("Invalid bounds ($l, $u) for row $i")
            end
        end

        # @info "Non-zeros after adding slacks" nz aI aJ aV
        
        # Populate objective, coefficients and bounds
        nz = 0
        nvarupb = 0
        nfree = 0
        for (j, (l, u)) in enumerate(zip(pb.lvar, pb.uvar))
            col = pb.acols[j]
            if l == -Inf && u == Inf
                # Free variable

                # Positive part
                c[j + nfree] = pb.obj[j]
                for (i, v) in zip(col.nzind, col.nzval)
                    nz += 1
                    aI[nz] = i
                    aJ[nz] = j + nfree
                    aV[nz] = v
                end

                # Negative part
                c[j + nfree + 1] = -pb.obj[j]
                for (i, v) in zip(col.nzind, col.nzval)
                    nz += 1
                    aI[nz] = i
                    aJ[nz] = j + nfree + 1
                    aV[nz] = -v
                end
                nfree += 1

            elseif l == -Inf && isfinite(u)
                # `xj ⩽ uj` is flipped to `xj = uj - xj_, xj_ ⩾ 0`
                # Thus, `aij xj == b` becomes `-aij*xj_ == b - aij*u`
                
                # Flip column and push new lower bound to zero
                for (i, v) in zip(col.nzind, col.nzval)
                    nz += 1
                    aI[nz] = i
                    aJ[nz] = j + nfree
                    aV[nz] = -v

                    b[i] -= v * u
                end

                # Populate objective
                # `cj xj` becomes `-cj xj_ + cj uj`
                c[j + nfree] = - pb.obj[j]
                c0 += pb.obj[j] * u

            elseif isfinite(l) && u == Inf
                # `xj ⩾ lj` becomes `xj = lj + xj_, xj_ ⩾ 0`

                # Push lower bound to zero
                for (i, v) in zip(col.nzind, col.nzval)
                    nz += 1
                    aI[nz] = i
                    aJ[nz] = j + nfree
                    aV[nz] = v

                    b[i] -= l * v
                end

                # Populate objective
                c[j + nfree] = pb.obj[j]
                c0 += pb.obj[j] * l

            elseif isfinite(l) && isfinite(u)
                # `l ⩽ x ⩽ u` becomes `0 ⩽ x_ ⩽ u-l`
                nvarupb += 1

                # Push lower bound to zero
                for (i, v) in zip(col.nzind, col.nzval)
                    nz += 1
                    aI[nz] = i
                    aJ[nz] = j + nfree
                    aV[nz] = v

                    b[i] -= l * v
                end

                # Populate objective
                c[j + nfree] = pb.obj[j]
                c0 += pb.obj[j] * l

                # Record upper bound
                uind[nvarupb] = j + nfree
                uval[nvarupb] = u - l
            else
                error("Invalid bounds ($l, $u) for variable $j")
            end
        end
        # @info "Non-zeros while populating matrix" nz
        # @info "A" aI aJ aV
        # If problem is maximization, flip objective
        if !pb.objsense
            c .= -c
            c0 = -c0
        end

        # Update dimensions
        nvar += nfree + nslack
        nupb = nvarupb + nslackupb

        # Build matrix
        A = construct_matrix(params.MatrixOptions.Ta, ncon, nvar, aI, aJ, aV; params.MatrixOptions.options...)

        # TODO: setup linear solver here
        
        return HSDSolver{Tv}(params, ncon, nvar, nupb, A, b, pb.objsense, c, c0, uind, uval)
    end

end

include("./hsd_step.jl")


"""
    compute_residuals!(::HSDSolver, res, pt, A, b, c, uind, uval)

In-place computation of primal-dual residuals at point `pt`.
"""
# TODO: check whether having just hsd as argument makes things slower
# TODO: Update solution status
function compute_residuals!(hsd::HSDSolver{T, Tv}, dat::IPMData{T, Tv, Tb, Ta}
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}}

    pt, res = hsd.pt, hsd.res

    # Primal residual
    # rp = t*b - A*x
    axpby!(pt.τ, dat.b, zero(T), res.rp)
    mul!(res.rp, dat.A, pt.x, -one(T), one(T))

    # Lower-bound residual
    # rl_j = τ*l_j - (x_j - xl_j)  if l_j ∈ R
    #      = 0                     if l_j = -∞
    @. res.rl = - pt.x + pt.xl + pt.τ * dat.l
    res.rl .*= dat.lflag

    # Upper-bound residual
    # ru_j = τ*u_j - (x_j + xu_j)  if u_j ∈ R
    #      = 0                     if u_j = +∞
    @. res.ru = - pt.x - pt.xu + pt.τ * dat.u
    res.ru .*= dat.uflag

    # Dual residual
    # rd = t*c - A'y - zl + zu
    axpby!(pt.τ, dat.c, zero(T), res.rd)
    mul!(res.rd, transpose(dat.A), pt.y, -one(T), one(T))
    @. res.rd += pt.zu .* dat.uflag - pt.zl .* dat.lflag

    # Gap residual
    # rg = c'x - (b'y + l'zl - u'zu) + k
    res.rg = pt.κ + (dot(dat.c, pt.x) - (
        dot(dat.b, pt.y)
        + dot(dat.l .* dat.lflag, pt.zl)
        - dot(dat.u .* dat.uflag, pt.zu)
    ))

    # Residuals norm
    res.rp_nrm = norm(res.rp, Inf)
    res.rl_nrm = norm(res.ru, Inf)
    res.ru_nrm = norm(res.ru, Inf)
    res.rd_nrm = norm(res.rd, Inf)
    res.rg_nrm = norm(res.rg, Inf)

    # Compute primal and dual bounds
    hsd.primal_bound_unscaled = dot(dat.c, pt.x) + pt.τ * dat.c0
    hsd.primal_bound_scaled   = hsd.primal_bound_unscaled / pt.τ
    hsd.dual_bound_unscaled   = pt.τ * dat.c0 + (
        dot(dat.b, pt.y)
        + dot(dat.l .* dat.lflag, pt.zl)
        - dot(dat.u .* dat.uflag, pt.zu)
    )
    hsd.dual_bound_scaled     = hsd.dual_bound_unscaled / pt.τ

    return nothing
end


"""
    update_solver_status!()

Update status and return true if solver should stop.
"""
function update_solver_status!(
    hsd::HSDSolver{T, Tv}, dat::IPMData{T, Tv, Tb, Ta},
    ϵp::T, ϵd::T, ϵg::T, ϵi::T
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}}
    hsd.solver_status = Trm_Unknown

    pt, res = hsd.pt, hsd.res

    ρp = max(
        res.rp_nrm / (pt.τ * (one(T) + norm(dat.b, Inf))),
        res.rl_nrm / (pt.τ * (one(T) + norm(dat.l .* dat.lflag, Inf))),
        res.ru_nrm / (pt.τ * (one(T) + norm(dat.u .* dat.uflag, Inf)))
    )
    ρd = res.rd_nrm / (pt.τ * (one(T) + norm(dat.c, Inf)))
    ρg = abs(hsd.primal_bound_unscaled - hsd.dual_bound_unscaled) / (pt.τ + abs(hsd.dual_bound_unscaled))

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
    
    # Check for infeasibility certificates
    if max(
        norm(dat.A * pt.x, Inf),
        norm((pt.x - pt.xl) .* dat.lflag, Inf),
        norm((pt.x + pt.xu) .* dat.uflag, Inf)
    ) * (norm(dat.c, Inf) / max(1, norm(dat.b, Inf))) < - ϵi * dot(dat.c, pt.x)
        # Dual infeasible, i.e., primal unbounded
        hsd.primal_status = Sln_InfeasibilityCertificate
        hsd.solver_status = Trm_DualInfeasible
        return nothing
    end

    δ = dat.A' * pt.y + (pt.zl .* dat.lflag) - (pt.zu .* dat.uflag)
    if norm(δ, Inf) * norm(dat.b, Inf) / (max(1, norm(dat.c, Inf)))  < (dot(dat.b, pt.y) + dot(dat.l .* dat.lflag, pt.zl)- dot(dat.u .* dat.uflag, pt.zu)) * ϵi
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
function ipm_optimize!(hsd::HSDSolver{T, Tv, Tk},
    dat::IPMData{T, Tv, Tb, Ta},
    params::Parameters{T}
) where{T, Tv<:AbstractVector{T}, Tb<:AbstractVector{Bool}, Ta<:AbstractMatrix{T}, Tk<:AbstractKKTSolver{T}}

    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function

    # Initialization
    TimerOutputs.reset_timer!(hsd.timer)
    tstart = time()
    hsd.niter = 0

    # Print information about the problem
    if params.OutputLevel > 0
        @printf "\nOptimizer info\n"
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

    # TODO: set starting point
    # Q: should we allocate memory for `pt` here?
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
        @timeit hsd.timer "Residuals" compute_residuals!(hsd, dat)

        update_mu!(hsd.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if params.OutputLevel > 0
            # Display log
            @printf "%4d" hsd.niter
            
            # Objectives
            ϵ = dat.objsense ? one(T) : -one(T)
            @printf "  %+14.7e" ϵ * hsd.primal_bound_scaled
            @printf "  %+14.7e" ϵ * hsd.dual_bound_scaled
            
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
        @timeit hsd.timer "update status" update_solver_status!(hsd, dat,
            params.BarrierTolerancePFeas,
            params.BarrierToleranceDFeas,
            params.BarrierToleranceRGap,
            params.BarrierToleranceIFeas
        )

        if (
            hsd.solver_status == Trm_Optimal
            || hsd.solver_status == Trm_PrimalInfeasible
            || hsd.solver_status == Trm_DualInfeasible
        )
            break
        elseif hsd.niter >= params.BarrierIterationsLimit 
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
            @timeit hsd.timer "Step" compute_step!(hsd, dat, params)
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