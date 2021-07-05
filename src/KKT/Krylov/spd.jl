"""
    SPDSolver
"""
mutable struct SPDSolver{T,V,Ta,KL,KS} <: AbstractKrylovSolver{T}
    # Problem data
    m::Int
    n::Int
    A::Ta

    # IPM-related workspace
    θ::Vector{T}
    regP::Vector{T}
    regD::Vector{T}

    # Krylov-related workspace
    D::Diagonal{T,V}
    Rd::Diagonal{T,V}
    ξ::V
    opK::KL

    # Krylov solver & related options
    atol::T
    rtol::T
    krylov_solver::KS

    # TODO: preconditioner
end

backend(kkt::SPDSolver) = "$(typeof(kkt.krylov_solver))"
linear_system(kkt::SPDSolver) = "K1"

function setup(A, ::K1, backend::Backend{KS,V}) where{KS<:Union{_KRYLOV_SPD,_KRYLOV_SID},V}
    Ta = typeof(A)
    T = eltype(A)
    T == eltype(V) || error("eltype(A)=$T incompatible with eltype of Krylov vector storage $V.")

    m, n = size(A)

    # Workspace
    θ = ones(T, n)
    regP = ones(T, n)
    regD = ones(T, m)

    D  = Diagonal(V(undef, n))
    Rd = Diagonal(V(undef, m))
    ξ  = V(undef, m)

    # Define linear operator for normal equations system
    # We need to allocate one temporary vector
    # This linear operator is symmetric definite positive,
    #   so we only need to define
    #   u ⟵ α (A×D×A'+Rd) × w + β u
    #   i.e., u = α(ADA')×w + (αRd)×w + βu
    vtmp = V(undef, n)
    opK = LO.LinearOperator(T, m, m, true, true,
        (u, w, α, β) -> begin
            mul!(vtmp, A', w)
            lmul!(D, vtmp)
            mul!(u, A, vtmp, α, β)
            mul!(u, Rd, w, α, one(T))
            u
        end
    )

    # Allocate Krylov solver's workspace
    atol = sqrt(eps(T))
    rtol = sqrt(eps(T))
    krylov_solver = KS(m, m, V)

    return SPDSolver{T,V,Ta,typeof(opK),typeof(krylov_solver)}(
        m, n, A,
        θ, regP, regD,
        D, Rd, ξ,
        opK,
        atol, rtol,
        krylov_solver
    )
end

function update!(kkt::SPDSolver, θ, regP, regD)

    copyto!(kkt.θ, θ)
    copyto!(kkt.regP, regP)
    copyto!(kkt.regD, regD)

    copyto!(kkt.Rd.diag, regD)
    copyto!(kkt.D.diag, inv.(kkt.θ .+ kkt.regP))

    return nothing
end

function solve!(dx, dy, kkt::SPDSolver{T}, ξp, ξd) where{T}
    m, n = kkt.m, kkt.n

    copyto!(kkt.ξ, ξp)

    mul!(kkt.ξ, kkt.A, kkt.D * ξd, true, true)

    # Solve the normal equations
    _krylov!(kkt.krylov_solver, kkt.opK, kkt.ξ; atol=kkt.atol, rtol=kkt.rtol)
    copyto!(dy, kkt.krylov_solver.x)

    # Recover dx
    copyto!(dx, ξd)
    mul!(dx, kkt.A', dy, one(T), -one(T))
    lmul!(kkt.D, dx)

    # TODO: iterative refinement (?)
    return nothing
end
