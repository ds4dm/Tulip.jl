"""
    test_linalg(A, b, c, uind, uval, h)

Verify that the linear algebra operations are properly defined for the input
    data structures.
"""
function test_linalg(A::AbstractMatrix{Tv}) where{Tv<:Real}

    # Dimension check
    m = size(A, 1)
    n = size(A, 2)

    # Matrix-vector multiplication
    x = zeros(Tv, n)
    y = zeros(Tv, m)

    # matrix-vector multiplication
    A * x;
    transpose(A) * y;
    mul!(y, A, x)

    # Cholesky factorization
    θ = Tv(2) .* ones(Tv, n)

    ls = TLP.TLPLinearSolver(A)
    TLP.TLPLinearAlgebra.update_linear_solver(ls, θ)

    # solve linear system
    dx = zeros(Tv, n)
    dy = zeros(Tv, m)
    ξp = ones(Tv, m)
    ξd = ones(Tv, n)
    TLP.TLPLinearAlgebra.solve_augmented_system!(
        dx, dy, ls, A, θ, ξp, ξd
    )

    # TODO: check that solution is approximately correct
    rp = norm(A  * dx - ξp, Inf)
    rd = norm(A' * dy - dx ./ θ - ξd, Inf)

    @test rp <= sqrt(eps(Tv))
    @test rd <= sqrt(eps(Tv))

    return true
end

@testset "LinearAlgebra" begin

    # Test specific data structures
    include("dense.jl")     # Dense matrices
    include("sparse.jl")    # SparseMatrixCSC
    
    # include("unitBlockAngular.jl") # Specialized unit block-angular

end