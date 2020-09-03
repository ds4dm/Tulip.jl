import Krylov
const LO = Krylov.LinearOperators

abstract type KrylovSolver{T} <: AbstractKKTSolver{T} end

# ==============================================================================
#   KrylovPDSolver: 
# ==============================================================================

"""
    KrylovPDSolver{T, F, Tv, Ta}

Solves normal equations system with an iterative method `f::F`.

```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(
    KrylovPDSolver, method=Krylov.cg
)
```
"""
mutable struct KrylovPDSolver{T, F, Tv, Ta} <: KrylovSolver{T}
    m::Int
    n::Int

    f::F

    # Problem data
    A::Ta     # Constraint matrix
    θ::Tv     # scaling
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization

    function KrylovPDSolver(f::Function, A::Ta) where{Ta <: AbstractMatrix}
        F = typeof(f)
        m, n = size(A)
        T = eltype(A)

        return new{T, F, Vector{T}, Ta}(
            m, n,
            f,
            A, ones(T, n), ones(T, n), ones(T, m)
        )
    end
end

setup(::Type{KrylovPDSolver}, A; method) = KrylovPDSolver(method, A)

backend(kkt::KrylovPDSolver) = "Krylov ($(kkt.f))"
linear_system(::KrylovPDSolver) = "Normal equations"


"""
    update!(kkt::KrylovPDSolver, θ, regP, regD)

Update diagonal scaling ``θ`` and primal-dual regularizations.
"""
function update!(
    kkt::KrylovPDSolver{Tv},
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
    solve!(dx, dy, kkt::KrylovPDSolver, ξp, ξd)

Solve the normal equations system using an iterative method.
"""
function solve!(
    dx::Vector{Tv}, dy::Vector{Tv},
    kkt::KrylovPDSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = kkt.m, kkt.n

    # Setup
    d = one(Tv) ./ (kkt.θ .+ kkt.regP)
    D = Diagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ kkt.A * (D * ξd)

    # Form linear operator S = A * D *A' + Rd
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
    _dy, stats = kkt.f(opS, ξ_)
    dy .= _dy

    # Recover dx
    dx .= D * (kkt.A' * dy - ξd)

    return nothing
end


# ==============================================================================
#   KrylovSQDSolver: 
# ==============================================================================

"""
    KrylovSQDSolver{T, F, Tv, Ta}

Solves the augmented system with an iterative method `f::F`.

```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(
    KrylovSQDSolver, method=Krylov.trimr
)
```
"""
mutable struct KrylovSQDSolver{T, F, Tv, Ta} <: KrylovSolver{T}
    m::Int
    n::Int

    f::F

    # Problem data
    A::Ta     # Constraint matrix
    θ::Tv     # scaling
    regP::Tv  # primal regularization
    regD::Tv  # dual regularization

    function KrylovSQDSolver(f::Function, A::Ta) where{Ta <: AbstractMatrix}
        F = typeof(f)
        m, n = size(A)
        T = eltype(A)
        return new{T, F, Vector{T}, Ta}(
            m, n,
            f,
            A, ones(T, n), ones(T, n), ones(T, m)
        )
    end
end

setup(::Type{KrylovSQDSolver}, A; method=Krylov.trimr) = KrylovSQDSolver(method, A)

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

Solve the normal equations system using an iterative method.
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
    _dy, _dx, stats = kkt.f(kkt.A, ξp, ξd, N=D, τ = kkt.regD[1])

    # @info stats.residuals
    dx .= _dx
    dy .= _dy

    return nothing
end