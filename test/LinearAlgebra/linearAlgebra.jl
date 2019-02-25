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
    F = Tulip.TLPLinearAlgebra.factor_normaleq(A, d)  # in-place update

    # solve linear system
    y = F \ b

    return true
end



# SparseCSC matrix
A = hcat(sparse(1.0I, 2, 2), sparse(1.0I, 2, 2))
m = A.m
n = A.n
c = [1.1, 1.2, 1.3, 1.4]
b = [1.1, 0.9]
u = sparse([0.5, 0.4, 0.7, 1.1])
uind = u.nzind
uval = u.nzval
p = nnz(u)

@test test_linalg(
    A, b, c, uind, uval,
    zeros(n), zeros(p), zeros(m), zeros(n), zeros(p)
)

# BlockAngular matrix
nblocks = 2
cols = [ones(m, 1) for _ in 1:nblocks]
B = Matrix{Float64}(I, m, m)
A = Tulip.TLPLinearAlgebra.UnitBlockAngular(cols, B)
(m, n) = size(A)
c = [1.1, 1.2, 1.3, 1.4]
b = [1.1, 1.2, 1.3, 1.4]
u = sparse([1.0, 0.0, 2.0, 3.0])
uind = u.nzind
uval = u.nzval
p = nnz(u)

@test test_linalg(
    A, b, c, uind, uval,
    zeros(n), zeros(p), zeros(m), zeros(n), zeros(p)
)