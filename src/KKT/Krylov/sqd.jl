"""
    SQDSolver
"""
mutable struct SQDSolver{T,V,Ta,KS} <: AbstractKrylovSolver{T}
    # Problem data
    m::Int
    n::Int
    A::Ta

    # IPM-related workspace
    θ::Vector{T}
    regP::Vector{T}
    regD::Vector{T}

    # Krylov-related workspace
    Θp::Diagonal{T,V}
    Θp⁻¹::Diagonal{T,V}
    Θd::Diagonal{T,V}
    Θd⁻¹::Diagonal{T,V}
    ξp::V
    ξd::V

    # Krylov solver & related options
    atol::T
    rtol::T
    krylov_solver::KS

    # TODO: preconditioner
end

backend(kkt::SQDSolver) = "$(typeof(kkt.krylov_solver))"
linear_system(kkt::SQDSolver) = "K2"

function setup(A, ::K2, backend::Backend{KS,V}) where{KS<:KrylovSQD,V}
    Ta = typeof(A)
    T = eltype(A)
    T == eltype(V) || error("eltype(A)=$T incompatible with eltype of Krylov vector storage $V.")

    m, n = size(A)

    # Workspace
    θ = ones(T, n)
    regP = ones(T, n)
    regD = ones(T, m)

    Θp   = Diagonal(V(undef, n))
    Θp⁻¹ = Diagonal(V(undef, n))
    Θd   = Diagonal(V(undef, m))
    Θd⁻¹ = Diagonal(V(undef, m))
    ξp  = V(undef, m)
    ξd  = V(undef, n)

    # Allocate Krylov solver's workspace
    atol = sqrt(eps(T))
    rtol = sqrt(eps(T))
    krylov_solver = KS(m, n, V)

    return SQDSolver{T,V,Ta,typeof(krylov_solver)}(
        m, n, A,
        θ, regP, regD,
        Θp, Θp⁻¹, Θd, Θd⁻¹,
        ξp, ξd,
        atol, rtol,
        krylov_solver
    )
end

function update!(kkt::SQDSolver, θ, regP, regD)

    copyto!(kkt.θ, θ)
    copyto!(kkt.regP, regP)
    copyto!(kkt.regD, regD)

    copyto!(kkt.Θp.diag, -(kkt.θ .+ kkt.regP))
    copyto!(kkt.Θp⁻¹.diag, inv.(kkt.θ .+ kkt.regP))  # Θp⁻¹ will be negated by tricg/trimr
    copyto!(kkt.Θd.diag, kkt.regD)
    copyto!(kkt.Θd⁻¹.diag, inv.(kkt.regD))

    return nothing
end

function solve!(dx, dy, kkt::SQDSolver{T}, ξp, ξd) where{T}
    copyto!(kkt.ξp, ξp)
    copyto!(kkt.ξd, ξd)

    # Solve the augmented system
    _krylov!(kkt.krylov_solver, kkt.A, kkt.ξp, kkt.ξd;
        M=kkt.Θd⁻¹,
        N=kkt.Θp⁻¹,
        atol=kkt.atol,
        rtol=kkt.rtol
    )

    # Recover dx, dy
    copyto!(dx, kkt.krylov_solver.yₖ)
    copyto!(dy, kkt.krylov_solver.xₖ)

    # TODO: iterative refinement (?)
    return nothing
end
