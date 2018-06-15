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
    Base.LinAlg.At_mul_B(A, b)
    Base.LinAlg.A_mul_B!(y, A, x)

    # Cholesky factorization
    F = Tulip.symbolic_cholesky(A)
    d = 1.0 + rand(n)
    Tulip.Cholesky.cholesky!(A, d, F)  # in-place update

    # solve linear system
    y = F \ b

    return nothing
end



# SparseCSC matrix
srand(0)
m = 4
n = 16
c = rand(n)
b = rand(m)
u = sprand(n, 0.75)
uind = u.nzind
uval = u.nzval
p = nnz(u)
A = sprand(m, n, 1.0)
test_linalg(A, b, c, uind, uval, zeros(n), zeros(p), zeros(m), zeros(n), zeros(p))

# BlockAngular matrix
srand(0)
nblocks = 4
cols = [rand(m, 4) for _ in 1:nblocks]
A = Tulip.Cholesky.DenseBlockAngular(cols)
(m, n) = size(A)
c = rand(n)
b = rand(m)
u = sprand(n, 0.75)
uind = u.nzind
uval = u.nzval
p = nnz(u)

test_linalg(A, b, c, uind, uval, zeros(n), zeros(p), zeros(m), zeros(n), zeros(p))