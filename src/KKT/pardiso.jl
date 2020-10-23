import Pardiso

mutable struct MKLPardisoSQD <: AbstractKKTSolver{Float64}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # Problem data
    A::SparseMatrixCSC{Float64, Int}
    θ::Vector{Float64}
    regP::Vector{Float64}  # primal regularization
    regD::Vector{Float64}  # dual regularization

    # Left-hand side matrix
    S::SparseMatrixCSC{Float64, Int}

    # Linear solver
    ps::Pardiso.MKLPardisoSolver

    function MKLPardisoSQD(A::SparseMatrixCSC{Float64})

        m, n = size(A)
        θ = ones(n)

        # We store we lower-triangular of the matrix
        S = [
            spdiagm(0 => -θ)  spzeros(n, m);
            A spdiagm(0 => ones(m))
        ]

        ps = Pardiso.MKLPardisoSolver()

        # We use symmetric indefinite matrices
        Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)
        Pardiso.pardisoinit(ps)

        # Set number of threads
        Pardiso.set_nprocs!(ps, 1)

        # Do the analysis
        Pardiso.set_phase!(ps, Pardiso.ANALYSIS)
        Pardiso.pardiso(ps, S, ones(m+n))

        return new(m, n, A, θ, ones(Float64, n), ones(Float64, m), S, ps)

        return kkt
    end
end

mutable struct PardisoSQD <: AbstractKKTSolver{Float64}
    m::Int  # Number of rows
    n::Int  # Number of columns

    # Problem data
    A::SparseMatrixCSC{Float64, Int}
    θ::Vector{Float64}
    regP::Vector{Float64}  # primal regularization
    regD::Vector{Float64}  # dual regularization

    # Left-hand side matrix
    S::SparseMatrixCSC{Float64, Int}

    # Linear solver
    ps::Pardiso.PardisoSolver

    function PardisoSQD(A::SparseMatrixCSC{Float64})

        m, n = size(A)
        θ = ones(n)

        # We store we lower-triangular of the matrix
        S = [
            spdiagm(0 => -θ)  A';
            A spdiagm(0 => ones(m))
        ]

        ps = Pardiso.PardisoSolver()

        # We use symmetric indefinite matrices
        Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)
        Pardiso.pardisoinit(ps)

        S_pardiso = Pardiso.get_matrix(ps, S, :N)
        @assert istril(S_pardiso)

        # Use a direct method
        Pardiso.set_solver!(ps, 1)

        # Do the analysis
        Pardiso.set_phase!(ps, Pardiso.ANALYSIS)
        Pardiso.pardiso(ps, S_pardiso, ones(m+n))

        return new(m, n, A, θ, ones(Float64, n), ones(Float64, m), S_pardiso, ps)

        return kkt
    end
end

setup(::Type{MKLPardisoSQD}, A) = MKLPardisoSQD(A)
backend(::MKLPardisoSQD) = "MKLPardiso"
linear_system(::MKLPardisoSQD) = "Augmented system"

setup(::Type{PardisoSQD}, A) = PardisoSQD(A)
backend(::PardisoSQD) = "Pardiso"
linear_system(::PardisoSQD) = "Augmented system"

"""
    update!(kkt, θ, regP, regD)

Update LDLᵀ factorization of the augmented system.

Update diagonal scaling ``\\theta``, primal-dual regularizations, and re-compute
    the factorization.
Throws a `PosDefException` if matrix is not quasi-definite.
"""
function update!(
    kkt::Union{MKLPardisoSQD, PardisoSQD},
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

    m, n = kkt.m, kkt.n

    # Update diagonal scaling
    kkt.θ .= θ
    # Update regularizers
    kkt.regP .= regP
    kkt.regD .= regD

    # Update S.
    # S is stored as lower-triangular and only its diagonal changes.
    @inbounds for j in 1:kkt.n
        k = kkt.S.colptr[j]
        kkt.S.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.S.colptr[kkt.n+i]
        kkt.S.nzval[k] = regD[i]
    end

    # Compute numerical factorization
    Pardiso.set_phase!(kkt.ps, Pardiso.NUM_FACT)
    Pardiso.pardiso(kkt.ps, kkt.S, zeros(kkt.m + kkt.n))

    return nothing
end

"""
    solve!(dx, dy, kkt, ξp, ξd)

Solve the augmented system, overwriting `dx, dy` with the result.
"""
function solve!(
    dx::Vector{Float64}, dy::Vector{Float64},
    kkt::Union{MKLPardisoSQD, PardisoSQD},
    ξp::Vector{Float64}, ξd::Vector{Float64}
)
    m, n = kkt.m, kkt.n
    
    # Set-up right-hand side
    ξ = [ξd; ξp]

    # Solve augmented system
    d = zeros(m + n)
    Pardiso.set_phase!(kkt.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    Pardiso.pardiso(kkt.ps, d, kkt.S, ξ)

    # Recover dx, dy
    @views dx .= d[1:n]
    @views dy .= d[(n+1):(m+n)]

    return nothing
end