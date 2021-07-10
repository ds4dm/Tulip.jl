const CholmodSQD = CholmodSolver{Float64,K2}

linear_system(::CholmodSQD) = "Augmented system (K2)"

function setup(A::SparseMatrixCSC{Float64,Int}, ::K2, ::Backend)
    m, n = size(A)

    θ = ones(Float64, n)
    regP = ones(Float64, n)
    regD = ones(Float64, m)
    ξ = zeros(Float64, m+n)

    K = [
        spdiagm(0 => -θ)  A';
        spzeros(Float64, m, n) spdiagm(0 => ones(m))
    ]

    # TODO: Symbolic factorization only
    F = ldlt(Symmetric(K))

    return CholmodSolver{Float64,K2}(m, n, A, θ, regP, regD, K, F, ξ)
end

function update!(kkt::CholmodSQD, θ, regP, regD)
    m, n = kkt.m, kkt.n

    # Sanity checks
    length(θ)  == n || throw(DimensionMismatch(
        "length(θ)=$(length(θ)) but KKT solver has n=$n."
    ))
    length(regP) == n || throw(DimensionMismatch(
        "length(regP)=$(length(regP)) but KKT solver has n=$n"
    ))
    length(regD) == m || throw(DimensionMismatch(
        "length(regD)=$(length(regD)) but KKT solver has m=$m"
    ))

    copyto!(kkt.θ, θ)
    copyto!(kkt.regP, regP)
    copyto!(kkt.regD, regD)

    # Update KKT matrix
    # K is stored as upper-triangular, and only its diagonal is changed
    @inbounds for j in 1:kkt.n
        k = kkt.K.colptr[1+j] - 1
        kkt.K.nzval[k] = -kkt.θ[j] - regP[j]
    end
    @inbounds for i in 1:kkt.m
        k = kkt.K.colptr[1+kkt.n+i] - 1
        kkt.K.nzval[k] = regD[i]
    end

    ldlt!(kkt.F, Symmetric(kkt.K))
    return nothing
end

function solve!(dx, dy, kkt::CholmodSQD, ξp, ξd)
    m, n = kkt.m, kkt.n

    # Setup right-hand side
    @views copyto!(kkt.ξ[1:n], ξd)
    @views copyto!(kkt.ξ[(n+1):end], ξp)

    # Solve augmented system
    # CHOLMOD doesn't have in-place solve, so this line will allocate
    δ = kkt.F \ kkt.ξ

    # Recover dx, dy
    @views copyto!(dx, δ[1:n])
    @views copyto!(dy, δ[(n+1):end])

    # TODO: iterative refinement
    return nothing
end
