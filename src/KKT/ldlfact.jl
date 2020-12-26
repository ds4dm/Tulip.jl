import LDLFactorizations
const LDLF = LDLFactorizations

# ==============================================================================
#   LDLFactSQD
# ==============================================================================

"""
    LDLFactSQD{T}

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses LDLFactorizations.jl to compute an LDLᵀ factorization of the quasi-definite augmented system.

```julia
model.params.KKT.Factory = Tulip.Factory(LDLFactSQD)
```
"""
mutable struct LDLFactSQD{T<:Real} <: AbstractKKTSolver{T}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{T, Int}
    θ::Vector{T}
    regP::Vector{T}  # primal regularization
    regD::Vector{T}  # dual regularization

    # Left-hand side matrix
    S::SparseMatrixCSC{T, Int}

    # Factorization
    F::LDLF.LDLFactorization{T}

    # TODO: constructor with initial memory allocation
    function LDLFactSQD(A::AbstractMatrix{T}) where{T<:Real}
        m, n = size(A)
        θ = ones(T, n)

        S = [
            spdiagm(0 => -θ)  A';
            spzeros(T, m, n) spdiagm(0 => ones(T, m))
        ]

        # TODO: PSD-ness checks
        # TODO: symbolic factorization only
        F = LDLF.ldl_analyze(Symmetric(S))

        return new{T}(m, n, A, θ, ones(T, n), ones(T, m), S, F)
    end

end

setup(::Type{LDLFactSQD}, A) = LDLFactSQD(A)

backend(::LDLFactSQD) = "LDLFactorizations.jl"
linear_system(::LDLFactSQD) = "Augmented system"

"""
    update!(kkt, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update!(
    kkt::LDLFactSQD{T},
    θ::AbstractVector{T},
    regP::AbstractVector{T},
    regD::AbstractVector{T}
) where{T<:Real}
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

    # Update S.
    # S is stored as upper-triangular and only its diagonal changes.
    @inbounds for j in 1:kkt.n
        k = kkt.S.colptr[1+j] - 1
        kkt.S.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.S.colptr[1+kkt.n+i] - 1
        kkt.S.nzval[k] = regD[i]
    end

    # TODO: PSD-ness checks
    try
        LDLF.ldl_factorize!(Symmetric(kkt.S), kkt.F)
    catch err
        isa(err, LDLF.SQDException) && throw(PosDefException(-1))
        rethrow(err)
    end

    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve!(
    dx::Vector{T}, dy::Vector{T},
    kkt::LDLFactSQD{T},
    ξp::Vector{T}, ξd::Vector{T}
) where{T<:Real}
    m, n = kkt.m, kkt.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = kkt.F \ ξ

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # TODO: Iterative refinement
    return nothing
end
