using SparseArrays
using SuiteSparse.CHOLMOD

"""
    Cholmod <: LSBackend

Use CHOLMOD backend.

Options available:
* `Float64` only
* Augmented system with LDLᵀ factorization
* Normal equations with Cholesky factorization
"""
struct Cholmod <: LSBackend end

# ==============================================================================
#   SparseIndefLinearSolver
# ==============================================================================

"""
    SparseIndefLinearSolver

Linear solver for the 2x2 augmented system with ``A`` sparse.

Uses an LDLᵀ factorization of the quasi-definite augmented system.
"""
mutable struct SparseIndefLinearSolver <: AbstractLinearSolver{Float64}
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
    function SparseIndefLinearSolver(A::SparseMatrixCSC{Float64, Int})
        m, n = size(A)
        θ = ones(Float64, n)

        S = [
            spdiagm(0 => -θ)  A';
            spzeros(Float64, m, n) spdiagm(0 => ones(m))
        ]

        # Symbolic factorization
        # See code in src/cholmod.jl of SuiteSparse.jl for `ldlt` function
        F = ldlt(Symmetric(S))

        return new(m, n, A, θ, ones(Float64, n), ones(Float64, m), S, F)
    end

end

AbstractLinearSolver(
    ::Cholmod,
    ::AugmentedSystem,
    A::AbstractMatrix{Float64}
) = SparseIndefLinearSolver(sparse(A))

AbstractLinearSolver(
    ::Cholmod,
    ::DefaultSystem,
    A::AbstractMatrix{Float64}
) = SparseIndefLinearSolver(sparse(A))

backend(::SparseIndefLinearSolver) = "CHOLMOD"
linear_system(::SparseIndefLinearSolver) = "Augmented system"

"""
    update_linear_solver!(ls, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update_linear_solver!(
    ls::SparseIndefLinearSolver,
    θ::AbstractVector{Float64},
    regP::AbstractVector{Float64},
    regD::AbstractVector{Float64}
)
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

    # Update S. S is stored as upper-triangular and only its diagonal changes.
    @inbounds for j in 1:ls.n
        k = ls.S.colptr[1+j] - 1
        ls.S.nzval[k] = -ls.θ[j] - regP[j]
    end
    @inbounds for i in 1:ls.m
        k = ls.S.colptr[1+ls.n+i] - 1
        ls.S.nzval[k] = regD[i]
    end

    # `S` is first converted into CHOLMOD's internal data structure,
    # so this line allocates.
    ldlt!(ls.F, Symmetric(ls.S))

    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve_augmented_system!(
    dx::Vector{Float64}, dy::Vector{Float64},
    ls::SparseIndefLinearSolver,
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = ls.m, ls.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = ls.F \ ξ

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    # Iterative refinement
    # TODO:
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = ls.A*dx + ls.regD .* dy - ξp
    # resD = - dx .* (ls.θ + ls.regP) + ls.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    # ξ1 = [resD; resP]
    # d1 = ls.F \ ξ1

    # # Update search direction
    # @views dx .-= d1[1:n]
    # @views dy .-= d1[(n+1):(m+n)]

    return nothing
end

# ==============================================================================
#   SparsePosDefLinearSolver
# ==============================================================================

"""
    SparsePosDefLinearSolver

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
mutable struct SparsePosDefLinearSolver <: AbstractLinearSolver{Float64}
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
    function SparsePosDefLinearSolver(A::SparseMatrixCSC{Float64, Int})
        m, n = size(A)
        θ = ones(Float64, n)

        S = A * A' + spdiagm(0 => ones(m))

        # TODO: PSD-ness checks
        F = cholesky(Symmetric(S))

        return new(m, n, A, θ, zeros(Float64, n), ones(Float64, m), F)
    end

end

AbstractLinearSolver(
    ::Cholmod,
    ::NormalEquations,
    A::AbstractMatrix{Float64}
) = SparsePosDefLinearSolver(sparse(A))

backend(::SparsePosDefLinearSolver) = "CHOLMOD - Cholesky"
linear_system(::SparsePosDefLinearSolver) = "Normal equations"

"""
    update_linear_solver!(ls, θ, regP, regD)

Compute normal equation system matrix, and update the factorization.
"""
function update_linear_solver!(
    ls::SparsePosDefLinearSolver,
    θ::AbstractVector{Float64},
    regP::AbstractVector{Float64},
    regD::AbstractVector{Float64}
)
    
    # Sanity checks
    length(θ) == ls.n || throw(DimensionMismatch(
        "θ has length $(length(θ)) but linear solver is for n=$(ls.n)."
    ))
    length(regP) == ls.n || throw(DimensionMismatch(
        "regP has length $(length(regP)) but linear solver has n=$(ls.n)"
    ))
    length(regD) == ls.m || throw(DimensionMismatch(
        "regD has length $(length(regD)) but linear solver has m=$(ls.m)"
    ))

    ls.θ .= θ
    ls.regP .= regP  # Primal regularization is disabled for normal equations
    ls.regD .= regD

    # Re-compute factorization
    # D = (Θ^{-1} + Rp)^{-1}
    D = Diagonal(one(Float64) ./ (ls.θ .+ ls.regP))
    Rd = spdiagm(0 => ls.regD)
    S = ls.A * D * ls.A' + Rd

    # Update factorization
    cholesky!(ls.F, Symmetric(S), check=false)
    issuccess(ls.F) || throw(PosDefException(0))

    return nothing
end

"""
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve_augmented_system!(
    dx::Vector{Float64}, dy::Vector{Float64},
    ls::SparsePosDefLinearSolver,
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = ls.m, ls.n

    d = one(Float64) ./ (ls.θ .+ ls.regP)
    D = Diagonal(d)
    
    # Set-up right-hand side
    ξ_ = ξp .+ ls.A * (D * ξd)

    # Solve augmented system
    dy .= (ls.F \ ξ_)

    # Recover dx
    dx .= D * (ls.A' * dy - ξd)

    # TODO: Iterative refinement
    # * Max number of refine steps
    # * Check for residuals before refining
    # * Check whether residuals did improve
    # resP = ls.A * dx + ls.regD .* dy - ξp
    # resD = - D \ dx + ls.A' * dy - ξd
    # println("\n|resP| = $(norm(resP, Inf))\n|resD| = $(norm(resD, Inf))")

    
    return nothing
end
