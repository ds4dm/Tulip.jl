print("\tdenseBlockAngular.jl")

srand(0)
m = 4
n = 2
R = 8
u = [rand(m, n) for _ in 1:R]

# Constructor test
A = Tulip.Cholesky.DenseBlockAngular(u)
@test m == A.m
@test R*n == A.n
@test R == A.R
for r in 1:R
    @test u[r] == A.cols[:, A.colptr[r]:(A.colptr[r+1]-1)]
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


println("\tPassed.")