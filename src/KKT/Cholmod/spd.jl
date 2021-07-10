const CholmodSPD = CholmodSolver{Float64,K1}

linear_system(::CholmodSPD) = "Normal equations (K1)"

function setup(A::SparseMatrixCSC{Float64}, ::K1, ::Backend)
    m, n = size(A)

    θ = ones(Float64, n)
    regP = ones(Float64, n)
    regD = ones(Float64, m)
    ξ = zeros(Float64, m)

    # TODO: analyze + in-place A*D*A' product
    K = sparse(A * A') + spdiagm(0 => regD)

    # TODO: PSD-ness checks
    F = cholesky(Symmetric(K))

    return CholmodSolver{Float64,K1}(m, n, A, θ, regP, regD, K, F, ξ)
end

function update!(kkt::CholmodSPD, θ, regP, regD)
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

    # Form normal equations matrix
    # TODO: use in-place update of S
    D = inv(Diagonal(kkt.θ .+ kkt.regP))
    kkt.K = (kkt.A * D * kkt.A') + spdiagm(0 => kkt.regD)

    # Update factorization
    cholesky!(kkt.F, Symmetric(kkt.K), check=false)
    issuccess(kkt.F) || throw(PosDefException(0))

    return nothing
end

function solve!(dx, dy, kkt::CholmodSPD, ξp, ξd)
    m, n = kkt.m, kkt.n

    D = inv(Diagonal(kkt.θ .+ kkt.regP))
    copyto!(kkt.ξ, ξp)
    mul!(kkt.ξ, kkt.A, D * ξd, true, true)

    # Solve normal equations
    # CHOLMOD doesn't have in-place solve, so this line will allocate
    dy .= kkt.F \ kkt.ξ

    # Recover dx
    copyto!(dx, ξd)
    mul!(dx, kkt.A', dy, 1.0, -1.0)
    lmul!(D, dx)

    # TODO: iterative refinement
    return nothing
end
