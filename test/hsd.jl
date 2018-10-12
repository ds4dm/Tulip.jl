using SparseArrays
using LinearAlgebra
using Random

m, n = 2, 3
Random.seed!(0)
A = sparse(rand(m, n))

x = rand(n)
s = rand(n)
w = rand(1)
z = rand(1)
uind = [2]

θwz = w ./ z
θ = s ./ x
θ[uind] .+= θwz
θ .\= 1.0

F = cholesky(Symmetric(A*Diagonal(θ)*A'))

rp = rand(m)
rd = rand(n)
ru = rand(1)

dx = zeros(n)
dy = zeros(m)
dz = zeros(1)

Tulip.solve_augsys_hsd!(A, F, θ, θwz, uind, dx, dy, dz, rp, rd, ru)

@test norm(A*dx - rp, Inf) <= 1e-8
dz_ = zeros(n)
dz_[uind] .= dz
@test norm(-(s./x).*dx + transpose(A)*dy - dz_ - rd, Inf) <= 1e-8
@test norm(dx[uind] - θwz.\dz - ru, Inf) <= 1e-8