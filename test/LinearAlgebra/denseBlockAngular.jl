print("\tdenseBlockAngular.jl")
Base.BLAS.set_num_threads(1)

m = 2
n = 1
R = 2
u = [ones(m, n) for _ in 1:R]
B = eye(m)

# Constructors
A = Tulip.LinearAlgebra.DenseBlockAngular([], zeros(0, 0))
@test size(A) == (0, 0)

A = Tulip.LinearAlgebra.DenseBlockAngular(u)
@test A.m == m
@test A.n == R*n

A = Tulip.LinearAlgebra.DenseBlockAngular([], B)
@test size(A) == size(B)
@test A.m == m
@test A.n == m
@test length(A.blocks) == 0
@test A.colslink == B

A = Tulip.LinearAlgebra.DenseBlockAngular(u, B)
@test size(A) == (n*R+m, m+R)
@test A.m == m
@test A.n == n*R+m
@test A.blocks == u


# Base interface tests
A_ = [
    [1.0 0.0 0.0 0.0];
    [0.0 1.0 0.0 0.0];
    [1.0 1.0 1.0 0.0];
    [1.0 1.0 0.0 1.0]
]
for i in Base.OneTo(A.n)
    for j in Base.OneTo(A.m+R)
        @test A[i, j] == A_[i, j]
    end
end


# Matrix-Vector multiplication tests
x = ones(A.n)
y = A * x
Base.LinAlg.A_mul_B!(y, A, x)
y_ = A_ * x
@test y == y_
x = Base.LinAlg.At_mul_B(A, y)
x_ = Base.LinAlg.At_mul_B(A_, y_)
@test x == x_

# Factorization tests
F = Tulip.LinearAlgebra.cholesky(A, ones(A.n))
@test m == F.m
@test R == F.R
@test (n*R+m) == F.n
@test size(F) == size(A)
@test F.colptr[end] == (F.n+1)

# factor update
θ = rand(A.n)
Tulip.LinearAlgebra.cholesky!(A, θ, F)

# Left division tests
A_ = sparse(A)  # sparse representation of A
b = rand(m+R)

y = F \ b
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10

Base.LinAlg.A_ldiv_B!(y, F, b)
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10

y = copy(b)
Base.LinAlg.A_ldiv_B!(F, y)
err = maximum(abs.(A_ * (θ .* (A_' * y)) - b))
@test err < 10.0^-10

# # Tulip tests
# # create and solve model
# m, n, R = 2, 2, 4
# u = [1.0 - 2.0 * rand(m, n) for _ in 1:R]
# for r in 1:R
#     u[r][:, 1] = 0.0
# end
# B = eye(m)
# A = Tulip.LinearAlgebra.DenseBlockAngular(u, B)
# b = vcat(ones(R), zeros(m))
# c = rand(A.n)
# colub_ind = collect(1:(A.n))
# colub_val = 10.0 * ones(A.n)

# # solve model
# model = Tulip.Model(A, b, c, colub_ind, colub_val)
# model.env[:verbose] = 0
# Tulip.optimize!(model)

println("\tPassed.")