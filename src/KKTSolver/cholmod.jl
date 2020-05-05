using SparseArrays
using SuiteSparse.CHOLMOD

"""
    CholmodSolver

Uses CHOLMOD's factorization routines to solve the augmented system.

To use an LDLᵀ factorization of the augmented system
(see [`Cholmod_SymQuasDef`](@ref))
```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(CholmodSolver, normal_equations=false)
```

To use a Cholesky factorization of the normal equations system
(see [`Cholmod_SymPosDef`](@ref))
```julia
model.params.KKTOptions = Tulip.KKT.SolverOptions(CholmodSolver, normal_equations=true)
```

!!! warning
    CHOLMOD can only be used with `Float64` arithmetic.
"""
abstract type CholmodSolver <: AbstractKKTSolver{Float64} end

function CholmodSolver(
    A::AbstractMatrix{Float64};
    normal_equations::Bool=false,
    kwargs...
)
    if normal_equations
        return Cholmod_SymPosDef(A; kwargs...)
    else
        return Cholmod_SymQuasDef(A; kwargs...)
    end
end


# ==============================================================================
#   Cholmod_SymQuasDef
# ==============================================================================

"""
    Cholmod_SymQuasDef

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses an LDLᵀ factorization of the quasi-definite augmented system.
"""
mutable struct Cholmod_SymQuasDef <: CholmodSolver
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Float64}
    θ::Vector{Float64}
    regP::Vector{Float64}  # primal-dual regularization
    regD::Vector{Float64}  # dual regularization

    # Left-hand side matrix
    S::SparseMatrixCSC{Float64, Int}  # TODO: use `CHOLMOD.Sparse` instead

    # Factorization
    F::CHOLMOD.Factor{Float64}

    # TODO: constructor with initial memory allocation
    function Cholmod_SymQuasDef(A::AbstractMatrix{Float64})
        m, n = size(A)
        θ = ones(Float64, n)

        S = [
            spdiagm(0 => -θ)  A';
            spzeros(Float64, m, n) spdiagm(0 => ones(m))
        ]

        # TODO: Symbolic factorization only
        F = ldlt(Symmetric(S))

        return new(m, n, A, θ, ones(Float64, n), ones(Float64, m), S, F)
    end

end

backend(::Cholmod_SymQuasDef) = "CHOLMOD"
linear_system(::Cholmod_SymQuasDef) = "Augmented system"

"""
    update!(kkt, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update!(
    kkt::Cholmod_SymQuasDef,
    θ::AbstractVector{Float64},
    regP::AbstractVector{Float64},
    regD::AbstractVector{Float64}
)
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

    # Update S. S is stored as upper-triangular and only its diagonal changes.
    @inbounds for j in 1:kkt.n
        k = kkt.S.colptr[1+j] - 1
        kkt.S.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.S.colptr[1+kkt.n+i] - 1
        kkt.S.nzval[k] = regD[i]
    end

    # `S` is first converted into CHOLMOD's internal data structure,
    # so this line allocates.
    ldlt!(kkt.F, Symmetric(kkt.S))

    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve!(
    dx::Vector{Float64}, dy::Vector{Float64},
    kkt::Cholmod_SymQuasDef,
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = kkt.m, kkt.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = kkt.F \ ξ

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # Iterative refinement
    # TODO:
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = kkt.A*dx + kkt.regD .* dy - ξp
    # resD = - dx .* (kkt.θ + kkt.regP) + kkt.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    # ξ1 = [resD; resP]
    # d1 = kkt.F \ ξ1

    # # Update search direction
    # @views dx .-= d1[1:n]
    # @views dy .-= d1[(n+1):(m+n)]

    return nothing
end

# ==============================================================================
#   Cholmod_SymPosDef
# ==============================================================================

"""
    Cholmod_SymPosDef

Linear solver for the 2x2 augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [xi_d]
    [   A          Rd] [dy]   [xi_p]
```
with ``A`` sparse.

Uses a Cholesky factorization of the positive definite normal equations system
```
(A * ((Θ⁻¹ + Rp)⁻¹ * Aᵀ + Rd) dy = xi_p + A * (θ⁻¹ + Rp)⁻¹ * xi_d
                              dx = (Θ⁻¹ + Rp)⁻¹ * (Aᵀ * dy - xi_d)
```
"""
mutable struct Cholmod_SymPosDef <: CholmodSolver
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Float64, Int}
    θ::Vector{Float64}   # Diagonal scaling
    # Regularization
    regP::Vector{Float64}  # primal
    regD::Vector{Float64}  # dual

    # Factorization
    F::CHOLMOD.Factor{Float64}

    # Constructor and initial memory allocation
    # TODO: symbolic only + allocation
    function Cholmod_SymPosDef(A::AbstractMatrix{Float64})
        m, n = size(A)
        θ = ones(Float64, n)

        S = sparse(A * A') + spdiagm(0 => ones(m))

        # TODO: PSD-ness checks
        F = cholesky(Symmetric(S))

        return new(m, n, A, θ, zeros(Float64, n), ones(Float64, m), F)
    end

end

backend(::Cholmod_SymPosDef) = "CHOLMOD - Cholesky"
linear_system(::Cholmod_SymPosDef) = "Normal equations"

"""
    update!(kkt, θ, regP, regD)

Compute normal equation system matrix, and update the factorization.
"""
function update!(
    kkt::Cholmod_SymPosDef,
    θ::AbstractVector{Float64},
    regP::AbstractVector{Float64},
    regD::AbstractVector{Float64}
)
    
    # Sanity checks
    length(θ) == kkt.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(kkt.n)."
    ))
    length(regP) == kkt.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(kkt.n)"
    ))
    length(regD) == kkt.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(kkt.m)"
    ))

    kkt.θ .= θ
    kkt.regP .= regP  # Primal regularization is disabled for normal equations
    kkt.regD .= regD

    # Re-compute factorization
    # D = (Θ^{-1} + Rp)^{-1}
    D = Diagonal(one(Float64) ./ (kkt.θ .+ kkt.regP))
    Rd = spdiagm(0 => kkt.regD)
    S = kkt.A * D * kkt.A' + Rd

    # Update factorization
    cholesky!(kkt.F, Symmetric(S), check=false)
    issuccess(kkt.F) || throw(PosDefException(0))

    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve!(
    dx::Vector{Float64}, dy::Vector{Float64},
    kkt::Cholmod_SymPosDef,
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = kkt.m, kkt.n

    d = one(Float64) ./ (kkt.θ .+ kkt.regP)
    D = Diagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ kkt.A * (D * ξd)

    # Solve augmented system
    dy .= (kkt.F \ ξ_)

    # Recover dx
    dx .= D * (kkt.A' * dy - ξd)

    # TODO: Iterative refinement
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = kkt.A * dx + kkt.regD .* dy - ξp
    # resD = - D \ dx + kkt.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    
    return nothing
end
