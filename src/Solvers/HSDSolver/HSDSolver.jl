"""
    HSDSolver

Solver for the homogeneous self-dual algorithm.
"""
mutable struct HSDSolver{Tv<:Real} <: AbstractIPMSolver{Tv}

    # =================
    #   Problem data
    # =================
    pb::StandardForm{Tv}


    # ================
    #   Book-keping
    # ================
    niter::Int  # Number of IPM iterations
    solver_status::TerminationStatus  # Optimization status
    primal_status::SolutionStatus
    dual_status::SolutionStatus


    #=====================
        Working memory
    =====================#

    # ncon::Int    # Number of constraints
    # nvar::Int    # Number of variables
    # nvarub::Int  # Number of upper-bounded variables

    pt::Point{Tv}    # Current primal-dual iterate
    res::Residuals{Tv}  # Residuals at current iterate
    F::Factorization{Tv}

    # rxs::Tv  # rhs for complimentary products
    # rwz::Tv  # rhs for complimentary products
    # rtk::T          # rhs for homogeneous complimentary products

    # TODO: Constructor
end


include("./hsd_step.jl")


"""
    compute_residuals!(::HSDSolver, res, pt, A, b, c, uind, uval)

In-place computation of primal-dual residuals at point `pt`.
"""
# TODO: check whether having just hsd as argument makes things slower
function compute_residuals!(
    ::HSDSolver,
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

    return nothing
end


"""
    check_stopping_criterion()

Check stopping criteria and return `true` if solver should stop.
"""
function check_stopping_criterion(
    hsd::HSDSolver,
    res::Residuals, pt::Point,
    A, b, c, uind, uval,

)
    stop = true

    return stop
end


"""
    optimize!

"""
function optimize!(hsd::HSDSolver{Tv}, env::TulipEnv) where{Tv<:Real}

    # TODO: pre-check whether model needs to be re-optimized.
    # This should happen outside of this function

    # Initialization
    tstart = time()
    niter = 0

    # TODO: allocate space for factorization

    # IPM LOG
    if env.verbose.val != 0
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
        compute_residuals!(hsd,
            hsd.res, hsd.pt,
            hsd.pb.A, hsd.pb.b, hsd.pb.c, hsd.pb.uind, hsd.pb.uval
        )

        update_mu!(hsd.pt)

        # I.B - Log
        # TODO: Put this in a logging function
        ttot = time() - tstart
        if env.verbose.val !=0
            # Display log
            @printf "%4d" niter
            
            # Objectives
            @printf "  %+16.7e" dot(hsd.pb.c, hsd.pt.x) / hsd.pt.t
            @printf "%+16.7e" (dot(hsd.pb.b, hsd.pt.y) - dot(hsd.pb.uval, hsd.pt.z)) / hsd.pt.t
            
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
        niter < env.barrier_iter_max.val || break
        ttot < env.time_limit.val || break

        # TODO: step
        # For now, include the factorization in the step function
        # Q: should we use more arguments here?
        compute_step!(hsd, env)

        niter += 1

    end

    return nothing

end