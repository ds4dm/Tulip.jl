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
    d = 1.1 .* ones(n)
    F = Tulip.TLPLinearAlgebra.factor_normaleq(A, d)

    # solve linear system
    y = F \ b

    return true
end

@testset "LinearAlgebra" begin

    # Test specific data structures
    include("sparseMatrixCSC.jl")  # General sparse matrices (Julia native)
    include("unitBlockAngular.jl") # Specialized unit block-angular

end