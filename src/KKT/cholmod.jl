using SparseArrays
using SuiteSparse.CHOLMOD

"""
    CholmodSolver

Uses CHOLMOD's factorization routines to solve the augmented system.

To use an LDLᵀ factorization of the augmented system
(see [`CholmodSQD`](@ref))
```julia
model.params.KKT.Factory = Tulip.Factory(CholmodSolver, normal_equations=false)
```

To use a Cholesky factorization of the normal equations system
(see [`CholmodSPD`](@ref))
```julia
model.params.KKT.Factory = Tulip.Factory(CholmodSolver, normal_equations=true)
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
        return CholmodSPD(A; kwargs...)
    else
        return CholmodSQD(A; kwargs...)
    end
end


# ==============================================================================
#   CholmodSQD
# ==============================================================================

"""
    CholmodSQD

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses an LDLᵀ factorization of the quasi-definite augmented system.
"""
mutable struct CholmodSQD <: CholmodSolver
    m::Int  # Number of rows
    n::Int  # Number of columns

    # TODO: allow for user-provided ordering,
    #   and add flag (default to false) to know whether user ordering should be used

    # Problem data
    A::SparseMatrixCSC{Float64}
    θ::Vector{Float64}     # diagonal scaling
    regP::Vector{Float64}  # primal-dual regularization
    regD::Vector{Float64}  # dual regularization
    ξ::Vector{Float64}     # right-hand side of the augmented system

    # Left-hand side matrix
    S::SparseMatrixCSC{Float64, Int}  # TODO: use `CHOLMOD.Sparse` instead

    # Factorization
    F::CHOLMOD.Factor{Float64}

    # TODO: constructor with initial memory allocation
    function CholmodSQD(A::AbstractMatrix{Float64})
        m, n = size(A)
        θ = ones(Float64, n)
        regP = ones(Float64, n)
        regD = ones(Float64, m)
        ξ = zeros(Float64, n+m)

        S = [
            spdiagm(0 => -θ)  A';
            spzeros(Float64, m, n) spdiagm(0 => ones(m))
        ]

        # TODO: Symbolic factorization only
        F = ldlt(Symmetric(S))

        return new(m, n, A, θ, regP, regD, ξ, S, F)
    end

end

backend(::CholmodSQD) = "CHOLMOD"
linear_system(::CholmodSQD) = "Augmented system"

"""
    update!(kkt, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update!(
    kkt::CholmodSQD,
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
    kkt::CholmodSQD,
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = kkt.m, kkt.n
    
    # Set-up the right-hand side
    @inbounds for i in 1:n
        kkt.ξ[i] = ξd[i]
    end
    @inbounds for i in 1:m
        kkt.ξ[i+n] = ξp[i]
    end

    # Solve augmented system
    d = kkt.F \ kkt.ξ

    # Recover dx, dy
    @inbounds for i in 1:n
        dx[i] = d[i]
    end
    @inbounds for i in 1:m
        dy[i] = d[i+n]
    end

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
#   CholmodSPD
# ==============================================================================

"""
    CholmodSPD

Linear solver for the 2x2 augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
with ``A`` sparse.

Uses a Cholesky factorization of the positive definite normal equations system
```
(A * ((Θ⁻¹ + Rp)⁻¹ * Aᵀ + Rd) dy = ξp + A * (θ⁻¹ + Rp)⁻¹ * ξd
                              dx = (Θ⁻¹ + Rp)⁻¹ * (Aᵀ * dy - ξd)
```
"""
mutable struct CholmodSPD <: CholmodSolver
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
    function CholmodSPD(A::AbstractMatrix{Float64})
        m, n = size(A)
        θ = ones(Float64, n)

        S = sparse(A * A') + spdiagm(0 => ones(m))

        # TODO: PSD-ness checks
        F = cholesky(Symmetric(S))

        return new(m, n, A, θ, zeros(Float64, n), ones(Float64, m), F)
    end

end

backend(::CholmodSPD) = "CHOLMOD - Cholesky"
linear_system(::CholmodSPD) = "Normal equations"

"""
    update!(kkt, θ, regP, regD)

Compute normal equation system matrix, and update the factorization.
"""
function update!(
    kkt::CholmodSPD,
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
    kkt::CholmodSPD,
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
