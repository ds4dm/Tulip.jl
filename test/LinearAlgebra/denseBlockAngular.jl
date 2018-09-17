using Random
using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(1)

Random.seed!(0)

m = 1
n = 1
R = 2
u = [ones(m, n) for _ in 1:R]
B = Matrix{Float64}(I, m, m)

# Constructors
A = Tulip.TLPLinearAlgebra.DenseBlockAngular([], zeros(0, 0))
@test size(A) == (0, 0)

A = Tulip.TLPLinearAlgebra.DenseBlockAngular(u)
@test A.m == m
@test A.n == R*n

A = Tulip.TLPLinearAlgebra.DenseBlockAngular([], B)
@test size(A) == size(B)
@test A.m == m
@test A.n == m
@test length(A.blocks) == 0
@test A.colslink == B

A = Tulip.TLPLinearAlgebra.DenseBlockAngular(u, B)
@test size(A) == (n*R+m, m+R)
@test A.m == m
@test A.n == n*R+m
@test A.blocks == u


# Base interface tests
A_ = [
    [1.0 0.0 0.0];
    [0.0 1.0 0.0];
    [1.0 1.0 1.0]
]
for i in Base.OneTo(A.n)
    for j in Base.OneTo(A.m+R)
        @test A[i, j] == A_[i, j]
    end
end

# Matrix-Vector product
x = ones(size(A, 2))
y = zeros(size(A, 1))
x_ = copy(x)
y_ = copy(y)

# A * x
mul!(y, A, x)
y_ = A_ * x_
@test y == y_

# A' * y
y = ones(size(A, 1))
y_ = ones(size(A, 1))
x_ = A' * y_
mul!(x, transpose(A), y)
@test x == x_

# Factorization tests
F = Tulip.TLPLinearAlgebra.cholesky(A, ones(A.n))
@test m == F.m
@test R == F.R
@test (n*R+m) == F.n
@test size(F) == size(A)
@test F.colptr[end] == (F.n+1)

# factor update
θ = rand(A.n)
Tulip.TLPLinearAlgebra.cholesky!(A, θ, F)



# Left division tests
b = ones(m+R)

y = F \ b
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10

ldiv!(y, F, b)
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10

y = copy(b)
ldiv!(F, y)
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10