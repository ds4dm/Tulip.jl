"""
    HSDSolver

Solver for the homogeneous self-dual algorithm.
"""
mutable struct HSDSolver{Tv<:Real} <: AbstractIPMSolver{Tv}

    # =================
    #   Problem data
    # =================
    ncon::Int  # Number of linear constraints
    nvar::Int  # Number of variables
    nupb::Int  # Number of upper-bounded variables
    A::AbstractMatrix{Tv}  # Constraint matrix
    b::Vector{Tv}  # Right-hand side
    objsense::Bool  # true if min, false if max
    c::Vector{Tv}  # Objective coefficients
    c0::Tv  # Objective offset
    uind::Vector{Int}  # Indices of upper-bounded variables
    uval::Vector{Tv}   # Upper-bounds on variables

    # =================
    #   Book-keeping
    # =================
    niter::Int  # Number of IPM iterations
    solver_status::TerminationStatus  # Optimization status
    primal_status::SolutionStatus
    dual_status::SolutionStatus

    primal_bound_unscaled::Tv  # Unscaled primal bound c'x
    primal_bound_scaled::Tv  # Scaled primal bound (c'x) / t
    dual_bound_unscaled::Tv  # Unscaled dual bound b'y - u'z
    dual_bound_scaled::Tv  # Scaled dual bound (b'y - u'z) / t


    #=====================
        Working memory
    =====================#
    pt::Point{Tv}    # Current primal-dual iterate
    res::Residuals{Tv}  # Residuals at current iterate
    ls::AbstractLinearSolver{Tv}
    regP::Vector{Tv}  # primal regularization
    regD::Vector{Tv}  # dual regularization
    regG::Tv  # gap regularization

    # rxs::Tv  # rhs for complimentary products
    # rwz::Tv  # rhs for complimentary products
    # rtk::T          # rhs for homogeneous complimentary products

    # TODO: Constructor
    function HSDSolver{Tv}(
        params::Parameters{Tv},
        ncon::Int, nvar::Int, nupb::Int,
        A::AbstractMatrix{Tv}, b::Vector{Tv}, objsense::Bool, c::Vector{Tv}, c0::Tv,
        uind::Vector{Int}, uval::Vector{Tv}
    ) where{Tv<:Real}
        hsd = new{Tv}()

        hsd.ncon = ncon
        hsd.nvar = nvar
        hsd.nupb = nupb
        hsd.A = A
        hsd.b = b
        hsd.objsense = objsense
        hsd.c = c
        hsd.c0 = c0
        hsd.uind = uind
        hsd.uval = uval

        hsd.niter = 0
        hsd.solver_status = Trm_Unknown
        hsd.primal_status = SolutionStatus(0)
        hsd.dual_status = SolutionStatus(0)

        hsd.primal_bound_scaled = Tv(Inf)
        hsd.primal_bound_unscaled = Tv(Inf)
        hsd.dual_bound_scaled = Tv(-Inf)
        hsd.dual_bound_unscaled = Tv(-Inf)

        hsd.pt = Point{Tv}(ncon, nvar, nupb)
        hsd.res = Residuals(
            zeros(Tv, ncon), zeros(Tv, nupb),
            zeros(Tv, nvar), zero(Tv),
            zero(Tv), zero(Tv), zero(Tv), zero(Tv)
        )

        hsd.ls = AbstractLinearSolver(params.LinearSolverBackend, params.LinearSolverSystem, A)

        # Initial regularizations
        hsd.regP = ones(Tv, nvar)
        hsd.regD = ones(Tv, ncon)
        hsd.regG = one(Tv)

        return hsd
    end

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
        A = construct_matrix(params.MatrixType, ncon, nvar, aI, aJ, aV)

        # TODO: setup linear solver
        
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
function compute_residuals!(
    hsd::HSDSolver{Tv},
    res::Residuals{Tv}, pt::Point{Tv},
    A::AbstractMatrix{Tv}, b::Vector{Tv}, c::Vector{Tv}, c0::Tv,
    uind::Vector{Int}, uval::Vector{Tv}
) where{Tv<:Real}

    # Primal residual
    # ``rp = t*b - A*x``
    mul!(res.rp, A, pt.x)
    rmul!(res.rp, -oneunit(Tv))
    axpy!(pt.t, b, res.rp)

    # Upper-bound residual
    # ``ru = t*u - w - x``
    rmul!(res.ru, zero(Tv))
    axpy!(-oneunit(Tv), pt.w, res.ru)
    @views axpy!(-oneunit(Tv), pt.x[uind], res.ru)
    axpy!(pt.t, uval, res.ru)

    # Dual residual
    # ``rd = t*c - A'*y - s + z``
    mul!(res.rd, transpose(A), pt.y)
    rmul!(res.rd, -oneunit(Tv))
    axpy!(pt.t, c, res.rd)
    axpy!(-oneunit(Tv), pt.s, res.rd)
    @views axpy!(oneunit(Tv), pt.z, res.rd[uind])

    # Gap residual
    res.rg = dot(c, pt.x) - (dot(b, pt.y) - dot(uval, pt.z)) + pt.k

    # Residuals norm
    res.rp_nrm = norm(res.rp, Inf)
    res.ru_nrm = norm(res.ru, Inf)
    res.rd_nrm = norm(res.rd, Inf)
    res.rg_nrm = norm(res.rg, Inf)

    # Compute primal and dual bounds
    hsd.primal_bound_unscaled = dot(c, pt.x) + pt.t * c0
    hsd.primal_bound_scaled   = hsd.primal_bound_unscaled / pt.t
    hsd.dual_bound_unscaled   = dot(b, pt.y) - dot(uval, pt.z) + pt.t * c0
    hsd.dual_bound_scaled     = hsd.dual_bound_unscaled / pt.t

    return nothing
end


"""
    update_solver_status!()

Update status and return true if solver should stop.
"""
function update_solver_status!(
    hsd::HSDSolver{Tv},
    pt::Point{Tv}, res::Residuals{Tv},
    A, b, c, c0, uind, uval,
    ϵp, ϵd, ϵg, ϵi
) where{Tv<:Real}
    hsd.solver_status = Trm_Unknown

    ρp = max(
        res.rp_nrm / (pt.t * (oneunit(Tv) + norm(b, Inf))),
        res.ru_nrm / (pt.t * (oneunit(Tv) + norm(uval, Inf)))
    )
    ρd = res.rd_nrm / (pt.t * (oneunit(Tv) + norm(c, Inf)))
    ρg = abs(hsd.primal_bound_unscaled - hsd.dual_bound_unscaled) / (pt.t + abs(hsd.dual_bound_unscaled))

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
    if max(norm(A*pt.x, Inf), norm(pt.x[uind] + pt.w, Inf)) * (norm(c, Inf) / max(1, norm(b, Inf))) < - ϵi * dot(c, pt.x)
        # Dual infeasible, i.e., primal unbounded
        hsd.primal_status = Sln_InfeasibilityCertificate
        hsd.solver_status = Trm_DualInfeasible
        return nothing
    end

    δ = A'pt.y + pt.s
    δ[uind] .-= pt.z
    if norm(δ, Inf) * norm(b, Inf) / (max(1, norm(c, Inf)))  < (dot(b, pt.y) - dot(uval, pt.z)) * ϵi
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
function optimize!(hsd::HSDSolver{Tv}, params::Parameters{Tv}) where{Tv<:Real}

    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function

    # Initialization
    tstart = time()
    hsd.niter = 0

    # Print information about the problem
    if params.OutputLevel > 0
        @printf "\nOptimizer info\n"
        @printf "Linear solver options\n"
        @printf "  %-12s : %s\n" "Precision" "$Tv"
        @printf "  %-12s : %s\n" "Backend" KKT.backend(hsd.ls)
        @printf "  %-12s : %s\n" "System" KKT.linear_system(hsd.ls)
    end

    # IPM LOG
    if params.OutputLevel > 0
        @printf "\n%4s  %14s  %14s  %8s %8s %8s  %7s  %4s\n" "Itn" "PObj" "DObj" "PFeas" "DFeas" "GFeas" "Mu" "Time"
    end

    # TODO: set starting point
    # Q: should we allocate memory for `pt` here?
    hsd.pt.x  .= oneunit(Tv)
    hsd.pt.w  .= oneunit(Tv)
    hsd.pt.t   = oneunit(Tv)
    hsd.pt.y  .= zero(Tv)
    hsd.pt.s  .= oneunit(Tv)
    hsd.pt.z  .= oneunit(Tv)
    hsd.pt.k   = oneunit(Tv)
    hsd.pt.qp .= zero(Tv)
    hsd.pt.qd .= zero(Tv)
    update_mu!(hsd.pt)

    # Main loop
    # Iteration 0 corresponds to the starting point.
    # Therefore, there is no numerical factorization before the first log is printed.
    # If the maximum number of iterations is set to 0, the only computation that occurs
    # is computing the residuals at the initial point.
    while(true)

        # I.A - Compute residuals at current iterate
        compute_residuals!(
            hsd,
            hsd.res, hsd.pt,
            hsd.A, hsd.b, hsd.c, hsd.c0, hsd.uind, hsd.uval
        )

        update_mu!(hsd.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if params.OutputLevel > 0
            # Display log
            @printf "%4d" hsd.niter
            
            # Objectives
            ϵ = hsd.objsense ? one(Tv) : -one(Tv)
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
        update_solver_status!(
            hsd, hsd.pt, hsd.res,
            hsd.A, hsd.b, hsd.c, hsd.c0, hsd.uind, hsd.uval,
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
            compute_step!(hsd, params)
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