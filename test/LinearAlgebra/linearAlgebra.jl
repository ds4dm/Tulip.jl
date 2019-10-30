"""
    test_linalg(A, b, c, uind, uval, h)

Verify that the linear algebra operations are properly defined for the input
    data structures.
"""
function test_linalg(
    A::AbstractMatrix{Tv},
    b::AbstractVector{Tv},
    c::AbstractVector{Tv},
    uind::AbstractVector{Int},
    uval::AbstractVector{Tv},
    x::AbstractVector{Tv},
    w::AbstractVector{Tv},
    y::AbstractVector{Tv},
    s::AbstractVector{Tv},
    z::AbstractVector{Tv}
) where{Tv<:Real}

    # Dimension check
    m = size(A, 1)
    n = size(A, 2)
    n == size(c, 1) || throw(DimensionMismatch(""))
    m == size(b, 1) || throw(DimensionMismatch(""))
    n == size(x, 1) || throw(DimensionMismatch(""))
    m == size(y, 1) || throw(DimensionMismatch(""))
    n == size(s, 1) || throw(DimensionMismatch(""))
    size(uind, 1) == size(uval, 1) || throw(DimensionMismatch(""))

    # matrix-vector multiplication
    A * x;
    transpose(A) * y;
    mul!(y, A, x)

    # Cholesky factorization
    θ = 1.1 .* ones(n)

    ls = TLP.TLPLinearSolver(A)
    TLP.TLPLinearAlgebra.update_linear_solver(ls, θ)

    # solve linear system
    dx = zeros(n)
    dy = zeros(m)
    ξp = ones(m)
    ξd = ones(n)
    TLP.TLPLinearAlgebra.solve_augmented_system!(
        dx, dy, ls, A, θ, ξp, ξd
    )

    # TODO: check that solution is approximately correct
    rp = norm(A  * dx - ξp, Inf)
    rd = norm(A' * dy - dx ./ θ - ξd, Inf)

    @test rp <= 1e-8
    @test rd <= 1e-8

    return true
end

@testset "LinearAlgebra" begin

    # Test specific data structures
    include("sparse.jl")            # General sparse matrices (Julia native)
    include("dense.jl")             # Dense matrices
    # include("unitBlockAngular.jl") # Specialized unit block-angular

end