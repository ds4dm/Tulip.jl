import Krylov
const LO = Krylov.LinearOperators

"""
    AbstractKrylovSolver{T}

Abstract type for Kyrlov-based linear solvers.
"""
abstract type AbstractKrylovSolver{T} <: AbstractKKTSolver{T} end

# ==============================================================================
#   KrylovSPDSolver:
# ==============================================================================

"""
    KrylovSPDSolver{T, F, Tv, Ta}

Krylov-based **S**ymmetric **P**ositive-**D**efinite (SPD) linear solver.

Applies a Krylov method to the normal equations systems, then recovers a solution
to the augmented system.
The selected Krylov method must therefore handle symmetric positive-definite systems.
Suitable methods are CG or CR.

A `KrylovSPDSolver` is selected as follows
```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(
    KrylovSPDSolver, method=Krylov.cg,
    atol=1e-12, rtol=1e-12
)
```

The `method` argument is a function `f::F` whose signature must match
```julia
dy, _ = f(S, ξ; atol, rtol)
```
where `S`, `ξ` and `dy` are the normal equations system's
left- and right-hand side, and solution vector, respectively.
`S` may take the form of a matrix, or of a suitable linear operator.
"""
mutable struct KrylovSPDSolver{T, F, Tv, Ta} <: AbstractKrylovSolver{T}
    m::Int
    n::Int

    f::F     # Krylov function
    atol::T  # absolute tolerance
    rtol::T  # relative tolerance

    # Problem data
    A::Ta     # Constraint matrix
    θ::Tv     # scaling
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization

    function KrylovSPDSolver(f::Function, A::Ta;
        atol::T=sqrt(eps(T)),
        rtol::T=sqrt(eps(T))
    ) where{T, Ta <: AbstractMatrix{T}}
        F = typeof(f)
        m, n = size(A)

        return new{T, F, Vector{T}, Ta}(
            m, n,
            f, atol, rtol,
            A, ones(T, n), ones(T, n), ones(T, m)
        )
    end
end

setup(::Type{KrylovSPDSolver}, A; method, kwargs...) = KrylovSPDSolver(method, A; kwargs...)

backend(kkt::KrylovSPDSolver) = "Krylov ($(kkt.f))"
linear_system(::KrylovSPDSolver) = "Normal equations"


"""
    update!(kkt::KrylovSPDSolver, θ, regP, regD)

Update diagonal scaling ``θ`` and primal-dual regularizations.
"""
function update!(
    kkt::KrylovSPDSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:Real}
    # Sanity checks
    length(θ)  == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    # Update diagonal scaling
    kkt.θ .= θ
    # Update regularizers
    kkt.regP .= regP
    kkt.regD .= regD

    return nothing
end

"""
    solve!(dx, dy, kkt::KrylovSPDSolver, ξp, ξd)

Solve the normal equations system using the selected Krylov method.
"""
function solve!(
    dx::Vector{Tv}, dy::Vector{Tv},
    kkt::KrylovSPDSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = kkt.m, kkt.n

    # Setup
    d = one(Tv) ./ (kkt.θ .+ kkt.regP)
    D = Diagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ kkt.A * (D * ξd)

    # Form linear operator S = A * D * A' + Rd
    # Here we pre-allocate the intermediary and final vectors
    v1 = zeros(Tv, n)
    v2 = zeros(Tv, m)
    opS = LO.LinearOperator(Tv, m, m, true, true,
        w -> begin
            mul!(v1, kkt.A', w)
            v1 .*= d
            mul!(v2, kkt.A, v1)
            v2 .+= kkt.regD .* w
            v2
        end
    )

    # Solve normal equations
    _dy, stats = kkt.f(opS, ξ_, atol=kkt.atol, rtol=kkt.rtol)
    dy .= _dy

    # Recover dx
    dx .= D * (kkt.A' * dy - ξd)

    return nothing
end


# ==============================================================================
#   KrylovSIDSolver:
# ==============================================================================

"""
    KrylovSIDSolver{T, F, Tv, Ta}

Krylov-based **S**ymmetric **I**n**D**efinite (SID) linear solver.

Applies a Krylov method to solve the augmented system,
without exploiting its quasi-definiteness properties.
The selected Krylov method must therefore handle symmetric indefinite systems.
Suitable methods are MINRES or MINRES-QLP.

A `KrylovSIDSolver` is selected as follows
```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(
    KrylovSIDSolver, method=Krylov.minres,
    atol=1e-12, rtol=1e-12
)
```

The `method` argument is a function `f::F` whose signature must match
```julia
Δ, _ = f(K, ξ; atol, rtol)
```
where `K`, `ξ` and `Δ` are the augmented system's
left- and right-hand side, and solution vector, respectively.
`K` may take the form of a matrix, or of a suitable linear operator.
"""
mutable struct KrylovSIDSolver{T, F, Tv, Ta} <: AbstractKrylovSolver{T}
    m::Int
    n::Int

    f::F     # Krylov function
    atol::T  # absolute tolerance
    rtol::T  # relative tolerance

    # Problem data
    A::Ta     # Constraint matrix
    θ::Tv     # scaling
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization

    function KrylovSIDSolver(f::Function, A::Ta;
        atol::T=sqrt(eps(T)),
        rtol::T=sqrt(eps(T))
    ) where{T, Ta <: AbstractMatrix{T}}
        F = typeof(f)
        m, n = size(A)

        return new{T, F, Vector{T}, Ta}(
            m, n,
            f, atol, rtol,
            A, ones(T, n), ones(T, n), ones(T, m)
        )
    end
end

setup(::Type{KrylovSIDSolver}, A; method, kwargs...) = KrylovSIDSolver(method, A; kwargs...)

backend(kkt::KrylovSIDSolver) = "Krylov ($(kkt.f))"
linear_system(::KrylovSIDSolver) = "Augmented system"


"""
    update!(kkt::KrylovSIDSolver, θ, regP, regD)

Update diagonal scaling ``θ`` and primal-dual regularizations.
"""
function update!(
    kkt::KrylovSIDSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:Real}
    # Sanity checks
    length(θ)  == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    # Update diagonal scaling
    kkt.θ .= θ
    # Update regularizers
    kkt.regP .= regP
    kkt.regD .= regD

    return nothing
end

"""
    solve!(dx, dy, kkt::KrylovSIDSolver, ξp, ξd)

Solve the augmented system using the selected Krylov method.
"""
function solve!(
    dx::Vector{Tv}, dy::Vector{Tv},
    kkt::KrylovSIDSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = kkt.m, kkt.n

    # Setup
    A = kkt.A
    D = Diagonal(-(kkt.θ .+ kkt.regP))  # D = Θ⁻¹ + Rp

    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Form linear operator K = [-(Θ⁻¹ + Rp) Aᵀ; A Rd]
    # Here we pre-allocate the final vector
    z = zeros(Tv, m+n)
    opK = LO.LinearOperator(Tv, n+m, n+m, true, false,
        w -> begin
            @views z1 = z[1:n]
            @views z2 = z[n+1:n+m]

            @views w1 = w[1:n]
            @views w2 = w[n+1:n+m]

            # z1 = -(Θ⁻¹ + Rp) w1 + A'w2
            mul!(z1, D, w1, one(Tv), zero(Tv))  # z1 = -(Θ⁻¹ + Rp) w1
            mul!(z1, A', w2, one(Tv), one(Tv))  # z1 += A'w2

            # z2 = A w1 + Rd w2
            mul!(z2, A, w1, one(Tv), zero(Tv))   # z2 = A w1
            mul!(z2, Diagonal(kkt.regD), w2, one(Tv), one(Tv))  # z2 += Rd * w2

            z
        end
    )

    # Solve augmented system
    Δ, stats = kkt.f(opK, ξ, atol=kkt.atol, rtol=kkt.rtol)

    # Recover dx, dy
    dx .= Δ[1:n]
    dy .= Δ[n+1:n+m]

    return nothing
end


# ==============================================================================
#   KrylovSQDSolver: 
# ==============================================================================

"""
    KrylovSQDSolver{T, F, Tv, Ta}

Krylov-based **S**ymmetric **Q**uasi-**D**efinite (SQD) linear solver.

Applies a Krylov method to solve the augmented system, taking advantage of
its 2x2 block structure and quasi-definiteness.
The selected Krylov method must therefore handle 2x2 symmetric quasi-definite systems.
Suitable methods are TriCG and TriMR.

A `KrylovSIDSolver` is selected as follows
```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(
    KrylovSQDSolver, method=Krylov.trimr,
    atol=1e-12, rtol=1e-12
)
```

The `method` argument is a function `f::F` whose signature must match
```
dy, dx, _ = f(A, ξp, ξd; M, N, atol, rtol)
```
where the augmented system is of the form
```
[ -N⁻¹ Aᵀ ] [dx] = [ξd]
[  A   M⁻¹] [dy] = [ξp]
```
i.e., ``N = (Θ^{-1} + R_{p})^{-1}`` and ``M = R_{d}^{-1}``.
"""
mutable struct KrylovSQDSolver{T, F, Tv, Ta} <: AbstractKrylovSolver{T}
    m::Int
    n::Int
    
    f::F     # Krylov function
    atol::T  # absolute tolerance
    rtol::T  # relative tolerance

    # Problem data
    A::Ta     # Constraint matrix
    θ::Tv     # scaling
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization

    function KrylovSQDSolver(f::Function, A::Ta;
        atol::T=sqrt(eps(T)),
        rtol::T=sqrt(eps(T))
    ) where{T, Ta <: AbstractMatrix{T}}
        F = typeof(f)
        m, n = size(A)

        return new{T, F, Vector{T}, Ta}(
            m, n,
            f, atol, rtol,
            A, ones(T, n), ones(T, n), ones(T, m)
        )
    end
end

setup(::Type{KrylovSQDSolver}, A; method, kwargs...) = KrylovSQDSolver(method, A; kwargs...)

backend(kkt::KrylovSQDSolver) = "Krylov ($(kkt.f))"
linear_system(::KrylovSQDSolver) = "Augmented system"


"""
    update!(kkt::KrylovSQDSolver, θ, regP, regD)

Update diagonal scaling ``θ`` and primal-dual regularizations.
"""
function update!(
    kkt::KrylovSQDSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:Real}
    # Sanity checks
    length(θ)  == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    # Update diagonal scaling
    kkt.θ .= θ
    # Update regularizers
    kkt.regP .= regP
    kkt.regD .= regD

    return nothing
end

"""
    solve!(dx, dy, kkt::KrylovSQDSolver, ξp, ξd)

Solve the augmented system using the selected Kyrlov method.
"""
function solve!(
    dx::Vector{Tv}, dy::Vector{Tv},
    kkt::KrylovSQDSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = kkt.m, kkt.n

    # Setup
    d = one(Tv) ./ (kkt.θ .+ kkt.regP)
    D = Diagonal(d)

    # Solve augmented system
    # Currently assumes that Rd is a multiple of the identity matrix
    _dy, _dx, stats = kkt.f(kkt.A, ξp, ξd,
        N=D, τ = kkt.regD[1],
        atol=kkt.atol, rtol=kkt.rtol
    )

    dx .= _dx
    dy .= _dy

    return nothing
end