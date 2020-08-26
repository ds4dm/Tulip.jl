import Krylov
# import LinearOperators
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
    d = one(Float64) ./ (kkt.θ .+ kkt.regP)
    D = LO.opDiagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ kkt.A * (D * ξd)

    # Solve normal equations
    opA = LO.LinearOperator(kkt.A)
    opS = opA * D * opA' + LO.opDiagonal(kkt.regD)
    y, stats = kkt.f(opS, ξ_)

    # @info stats.residuals
    dy .= y

    # Recover dx
    dx .= D * (kkt.A' * dy - ξd)

    return nothing
end