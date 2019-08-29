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
    c::Vector{Tv}  # Objective coefficients
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
    F::Union{Nothing, Factorization{Tv}}

    # rxs::Tv  # rhs for complimentary products
    # rwz::Tv  # rhs for complimentary products
    # rtk::T          # rhs for homogeneous complimentary products

    # TODO: Constructor
    function HSDSolver{Tv}(ncon::Int, nvar::Int, nupb::Int,
        A::AbstractMatrix{Tv}, b::Vector{Tv}, c::Vector{Tv},
        uind::Vector{Int}, uval::Vector{Tv}
    ) where{Tv<:Real}
        hsd = new{Tv}()

        hsd.ncon = ncon
        hsd.nvar = nvar
        hsd.nupb = nupb
        hsd.A = A
        hsd.b = b
        hsd.c = c
        hsd.uind = uind
        hsd.uval = uval

        hsd.niter = 0
        hsd.solver_status = TerminationStatus(0)
        hsd.primal_status = SolutionStatus(0)
        hsd.dual_status = SolutionStatus(0)

        hsd.primal_bound_scaled = Tv(Inf)
        hsd.primal_bound_unscaled = Tv(Inf)
        hsd.dual_bound_scaled = Tv(-Inf)
        hsd.dual_bound_unscaled = Tv(-Inf)

        hsd.pt = Point(
            ncon, nvar, nupb,
            ones(Tv, nvar), ones(Tv, nupb), oneunit(Tv),
            ones(Tv, ncon), ones(Tv, nvar), ones(Tv, nupb), oneunit(Tv),
            oneunit(Tv)
        )
        hsd.res = Residuals(
            zeros(Tv, ncon), zeros(Tv, nupb),
            zeros(Tv, nvar), zero(Tv),
            zero(Tv), zero(Tv), zero(Tv), zero(Tv)
        )

        hsd.F = nothing

        return hsd
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
    A::AbstractMatrix{Tv}, b::Vector{Tv}, c::Vector{Tv},
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
    hsd.primal_bound_unscaled = dot(c, pt.x)
    hsd.primal_bound_scaled   = hsd.primal_bound_unscaled / pt.t
    hsd.dual_bound_unscaled   = dot(b, pt.y) - dot(uval, pt.z)
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
    A, b, c, uind, uval,
    ϵp, ϵd, ϵg, ϵi
) where{Tv<:Real}
    hsd.solver_status = TerminationStatus(0)
    
    # Check for optimal solution
    if (
        res.rp_nrm     <= ϵp * pt.t * (oneunit(Tv) + norm(b, Inf))
        && res.ru_nrm  <= ϵp * pt.t * (oneunit(Tv) + norm(uval, Inf))
        && res.rd_nrm  <=  ϵd * pt.t * (oneunit(Tv) + norm(c, Inf))
        && (abs(hsd.primal_bound_unscaled - hsd.dual_bound_unscaled) / (pt.t + abs(hsd.dual_bound_unscaled))) <= ϵg
    )
        # Optimality reached
        hsd.solver_status = TerminationStatus(1)
        return nothing
    end

    # Check for infeasibility certificates
    if (pt.μ < ϵi) && ((pt.t / pt.k) < ϵi)
        # Check for primal or dual infeasibility
        if dot(c, pt.x) < -ϵi
            # Dual infeasible
            hsd.solver_status = TerminationStatus(3)
            return nothing
        
        elseif (-dot(b, pt.y) - dot(uval, pt.z)) < -ϵi
            # Primal infeasible
            hsd.solver_status = TerminationStatus(2)
            return nothing
        else
            # Should never be reached
            error("Detected infeasibility, but")
        end
    end

    return nothing
end


"""
    optimize!

"""
function optimize!(hsd::HSDSolver{Tv}, env::Env{Tv}) where{Tv<:Real}

    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function

    # Initialization
    tstart = time()
    hsd.niter = 0

    # TODO: allocate space for factorization
    hsd.F = symbolic_cholesky(hsd.A)

    # IPM LOG
    if env.verbose != 0
        @printf "%4s  %16s%16s %9s%9s%9s  %7s  %4s\n" "Itn" "PObj" "DObj" "PFeas" "DFeas" "GFeas" "Mu" "Time"
        # println(" Itn    Primal Obj      Dual Obj        PFeas    DFeas    GFeas     Mu       Time")
    end


    # TODO: compute symbolic Cholesky
    # For this to be efficient, we need to now the type of A in the signature

    # TODO: set starting point
    # Q: should we allocate memory for `pt` here?
    hsd.pt.x .= oneunit(Tv)
    hsd.pt.w .= oneunit(Tv)
    hsd.pt.t  = oneunit(Tv)
    hsd.pt.y .= zero(Tv)
    hsd.pt.s .= oneunit(Tv)
    hsd.pt.z .= oneunit(Tv)
    hsd.pt.k  = oneunit(Tv)
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
            hsd.A, hsd.b, hsd.c, hsd.uind, hsd.uval
        )

        update_mu!(hsd.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if env.verbose !=0
            # Display log
            @printf "%4d" hsd.niter
            
            # Objectives
            @printf "  %+16.7e" dot(hsd.c, hsd.pt.x) / hsd.pt.t
            @printf "%+16.7e" (dot(hsd.b, hsd.pt.y) - dot(hsd.uval, hsd.pt.z)) / hsd.pt.t
            
            # Residuals
            @printf " %9.2e" max(hsd.res.rp_nrm, hsd.res.ru_nrm)
            @printf "%9.2e" hsd.res.rd_nrm
            @printf "%9.2e" hsd.res.rg_nrm

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
            hsd.A, hsd.b, hsd.c, hsd.uind, hsd.uval,
            env.barrier_tol_pfeas,
            env.barrier_tol_dfeas,
            env.barrier_tol_conv,
            env.barrier_tol_infeas
        )

        if (
            hsd.solver_status == TerminationStatus(1)     # Optimal
            || hsd.solver_status == TerminationStatus(2)  # Primal infeasible
            || hsd.solver_status == TerminationStatus(3)  # Dual infeasible
        )
            break
        elseif hsd.niter >= env.barrier_iter_max 
            hsd.solver_status = TerminationStatus(4)  # Iteration limit
            break
        elseif ttot >= env.time_limit
            hsd.solver_status = TerminationStatus(5)  # Iteration limit
            break
        end
        

        # TODO: step
        # For now, include the factorization in the step function
        # Q: should we use more arguments here?
        compute_step!(hsd, env)

        hsd.niter += 1

    end

    # TODO: print message based on termination status
    env.verbose == 1 && println("Solver exited with status $((hsd.solver_status))")

    return nothing

end