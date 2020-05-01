import LDLFactorizations
const LDLF = LDLFactorizations


"""
    LDLFact <: LSBackend

Use LDLFactorizations backend.

Options available:
* Any numerical type `T`
* Augmented system with LDLᵀ factorization
"""
struct LDLFact <: LSBackend end

# ==============================================================================
#   LDLFLinearSolver
# ==============================================================================

"""
    LDLFLinearSolver{Tv}

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses LDLFactorizations.jl to compute an LDLᵀ factorization of the quasi-definite augmented system.
"""
mutable struct LDLFLinearSolver{Tv<:Real} <: AbstractLinearSolver{Tv}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Tv, Int}
    θ::Vector{Tv}
    regP::Vector{Tv}  # primal regularization
    regD::Vector{Tv}  # dual regularization

    # Factorization
    F::LDLF.LDLFactorization{Tv}

    # TODO: constructor with initial memory allocation
    function LDLFLinearSolver(A::SparseMatrixCSC{Tv, Int}) where{Tv<:Real}
        m, n = size(A)
        θ = ones(Tv, n)

        S = [
            spdiagm(0 => -θ)  A';
            A spdiagm(0 => ones(m))
        ]

        # TODO: PSD-ness checks
        F = LDLF.ldl(S)
        return new{Tv}(m, n, A, θ, ones(Tv, n), ones(Tv, m), F)
    end

end

AbstractLinearSolver(
    ::LDLFact,
    ::AugmentedSystem,
    A::AbstractMatrix{Tv}
) where{Tv<:Real} = LDLFLinearSolver(sparse(A))

backend(::LDLFLinearSolver) = "LDLFactorizations.jl"
linear_system(::LDLFLinearSolver) = "Augmented system"

"""
    update_linear_solver!(ls, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update_linear_solver!(
    ls::LDLFLinearSolver{Tv},
    θ::AbstractVector{Tv},
    regP::AbstractVector{Tv},
    regD::AbstractVector{Tv}
) where{Tv<:Real}
    # Sanity checks
    length(θ)  == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(regP) == ls.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(ls.n)"
    ))
    length(regD) == ls.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(ls.m)"
    ))

    # Update diagonal scaling
    ls.θ .= θ
    # Update regularizers
    ls.regP .= regP
    ls.regD .= regD

    # Re-compute factorization
    # TODO: Keep S in memory, only change diagonal
    S = [
        spdiagm(0 => -ls.θ .- regP)  ls.A';
        ls.A spdiagm(0 => regD)
    ]

    # TODO: PSD-ness checks
    try
        ls.F = LDLF.ldl(S)
    catch err
        isa(err, LDLF.SQDException) && throw(PosDefException(-1))
        rethrow(err)
    end

    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve_augmented_system!(
    dx::Vector{Tv}, dy::Vector{Tv},
    ls::LDLFLinearSolver{Tv},
    ξp::Vector{Tv}, ξd::Vector{Tv}
) where{Tv<:Real}
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = ls.F \ ξ

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # TODO: Iterative refinement
    return nothing
end
