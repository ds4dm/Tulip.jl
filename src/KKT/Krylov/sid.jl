"""
    SIDSolver
"""
mutable struct SIDSolver{T,V,Ta,KL,KS} <: AbstractKrylovSolver{T}
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
    Θd::Diagonal{T,V}
    ξ::V
    opK::KL

    # Krylov solver & related options
    atol::T
    rtol::T
    krylov_solver::KS

    # TODO: preconditioner
end

backend(kkt::SIDSolver) = "$(typeof(kkt.krylov_solver))"
linear_system(kkt::SIDSolver) = "K2"

function setup(A, ::K2, backend::Backend{KS,V}) where{KS<:_KRYLOV_SID,V}
    Ta = typeof(A)
    T = eltype(A)
    T == eltype(V) || error("eltype(A)=$T incompatible with eltype of Krylov vector storage $V.")

    m, n = size(A)

    # Workspace
    θ = ones(T, n)
    regP = ones(T, n)
    regD = ones(T, m)

    Θp = Diagonal(V(undef, n))
    Θd = Diagonal(V(undef, m))
    ξ  = V(undef, m+n)

    # Define linear operator for the augmented system
    # This linear operator is symmetric indefinite
    opK = LO.LinearOperator(T, m+n, m+n, true, false,
        (u, w, α, β) -> begin
            @views u1 = u[1:n]
            @views u2 = u[(n+1):(m+n)]
            @views w1 = w[1:n]
            @views w2 = w[n+1:n+m]

            mul!(u1, Θp, w1, α, β)
            mul!(u1, A', w2, α, one(T))

            mul!(u2, A , w1, α, β)
            mul!(u2, Θd, w2, α, one(T))
            u
        end
    )

    # Allocate Krylov solver's workspace
    atol = sqrt(eps(T))
    rtol = sqrt(eps(T))
    krylov_solver = KS(m+n, m+n, V)

    return SIDSolver{T,V,Ta,typeof(opK),typeof(krylov_solver)}(
        m, n, A,
        θ, regP, regD,
        Θp, Θd, ξ,
        opK,
        atol, rtol,
        krylov_solver
    )
end

function update!(kkt::SIDSolver, θ, regP, regD)

    copyto!(kkt.θ, θ)
    copyto!(kkt.regP, regP)
    copyto!(kkt.regD, regD)

    copyto!(kkt.Θd.diag, regD)
    copyto!(kkt.Θp.diag, -(kkt.θ .+ kkt.regP))

    return nothing
end

function solve!(dx, dy, kkt::SIDSolver{T}, ξp, ξd) where{T}
    m, n = kkt.m, kkt.n

    @views copyto!(kkt.ξ[1:m], ξp)
    @views copyto!(kkt.ξ[(m+1):(m+n)], ξd)

    # Solve the augmented system
    _krylov!(kkt.krylov_solver, kkt.opK, kkt.ξ; atol=kkt.atol, rtol=kkt.rtol)

    # Recover dx, dy
    copyto!(dx, kkt.krylov_solver.x[1:n])
    copyto!(dy, kkt.krylov_solver.x[(n+1):(m+n)])

    # TODO: iterative refinement (?)
    return nothing
end
