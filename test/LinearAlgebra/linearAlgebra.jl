"""
    test_linalg(A, b, c, uind, uval, h)

Verifies that the linear algebra operations are properly defined for the matrix
    structure given in input.
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
    d = 1.0 .+ rand(n)
    F = Tulip.TLPLinearAlgebra.cholesky(A, d)  # in-place update

    # solve linear system
    y = F \ b

    return nothing
end



# SparseCSC matrix
using SparseArrays

Random.seed!(0)
A = hcat(sparse(1.0I, 2, 2), sparse(1.0I, 2, 2))
m = A.m
n = A.n
c = rand(n)
b = rand(m)
u = sprand(n, 1.0)
uind = u.nzind
uval = u.nzval
p = nnz(u)

test_linalg(A, b, c, uind, uval, zeros(n), zeros(p), zeros(m), zeros(n), zeros(p))

# BlockAngular matrix
Random.seed!(0)
nblocks = 2
cols = [ones(m, 1) for _ in 1:nblocks]
B = Matrix{Float64}(I, m, m)
A = Tulip.TLPLinearAlgebra.DenseBlockAngular(cols, B)
(m, n) = size(A)
c = rand(n)
b = rand(m)
u = sprand(n, 0.75)
uind = u.nzind
uval = u.nzval
p = nnz(u)

test_linalg(A, b, c, uind, uval, zeros(n), zeros(p), zeros(m), zeros(n), zeros(p))