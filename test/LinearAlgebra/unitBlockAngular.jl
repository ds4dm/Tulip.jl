using LinearAlgebra

m = 2
R = 2
B = [
    [1.0 1.0 1.0 1.0];
    [0.0 2.0 -3.0 2.0]
]
B0 = Matrix(1.0I, m, m)
blockidx = [1, 2, 1, 2]

d = [0.5, 0.2, 0.1, 9.9]
θ = [9.9, 0.9, 2.0, 3.0, 1.0, 5.0]

A = Tulip.TLPLinearAlgebra.UnitBlockAngular(B, R, blockidx, B0)
A_ = sparse(A)

x = ones(size(A, 2))
y = A*x
y_ = A_ * x

@test (y ≈ y_)

x = A' * y
x_ = A_' * y_

@test (x ≈ x_)


C = zeros(A.m, A.m)
B_ = hcat(A.B0, A.B)
C_ = A.B * Diagonal((d)) * A.B'

Tulip.TLPLinearAlgebra.rank_update!(0.0, C, 1.0, A.B, A.B_, d)
@test Symmetric(C) ≈ Symmetric(C_)

Tulip.TLPLinearAlgebra.rank_update!(0.0, C, 1.0, A.B, d)
@test Symmetric(C) ≈ Symmetric(C_)


F = Tulip.TLPLinearAlgebra.factor_normaleq(A, θ)
Tulip.TLPLinearAlgebra.factor_normaleq!(A, θ, F)

b = ones(A.M)
ldiv!(F, b)
r = ones(A.M) .- A* (θ .* (A'*b))  # residual

@test norm(r, Inf) <= 1e-10

A2 = Tulip.TLPLinearAlgebra.UnitBlockAngular([i*ones(2, 2) for i in 1:2], sparse(1.0I, 2, 2))
A3 = Tulip.TLPLinearAlgebra.UnitBlockAngular([i*ones(2, 2) for i in 1:2])