print("\tdenseBlockAngular.jl")
Base.BLAS.set_num_threads(1)

srand(0)
m = 24
n = 32
R = 16
u = [rand(m, n) for _ in 1:R]

# Constructor test
A = Tulip.LinearAlgebra.DenseBlockAngular(u)
@test m == A.m
@test R*n == A.n
@test R == A.R
for r in 1:R
    @test u[r] == A.cols[r]
end


# Base interface tests
@test size(A) == (m+R, n*R)
@test A[1, 1] == 1.0
@test A[1, end] == (R == 1 ? 1.0 : 0.0)
@test A[end, end] == u[end][end, end]
@test A[R+1, 1] == u[1][1, 1]


# Matrix-Vector multiplication tests
x = rand(A.n)
# these two should not throw any error
y = A * x
Base.LinAlg.A_mul_B!(y, A, x)
y_ = hcat(u...) * x
@test y[(R+1):end] == y_


# factorization tests
F = Tulip.LinearAlgebra.cholesky(A, ones(A.n))
@test m == F.m
@test R == F.R
@test n*R == F.n
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

# Tulip tests
# create and solve model
m, n, R = 2, 2, 4
u = [1.0 - 2.0 * rand(m, n) for _ in 1:R]
for r in 1:R
    u[r][:, 1] = 0.0
end
A = Tulip.LinearAlgebra.DenseBlockAngular(u)
b = vcat(ones(R), zeros(m))
c = rand(n*R)
colub_ind = collect(1:(n*R))
colub_val = 10.0 * ones(n*R)

# solve model
model = Tulip.Model(A, b, c, colub_ind, colub_val)
model.env.output_level = 0
Tulip.optimize!(model)


println("\tPassed.")