using SparseArrays
using LinearAlgebra
using Random

function test_augmented_system(A, F, θ, θxs, θwz, uind, ξp, ξd, ξu)

    dx = zeros(length(ξd))
    dy = zeros(length(ξp))
    dz = zeros(length(ξu))

    # Solve augmented system
    Tulip.solve_augsys_hsd!(A, F, θ, θwz, uind, dx, dy, dz, ξp, ξd, ξu)

    # Compute residuals and check numerical precision
    rp = A*dx .- ξp

    dz_ = zeros(length(ξd))
    dz_[uind] .= dz
    rd = - θxs .* dx .+ transpose(A)*dy .- dz_ .- ξd

    ru = dx[uind] .- θwz .\ dz - ξu

    @test norm(rp, Inf) <= 1e-12
    @test norm(rd, Inf) <= 1e-12
    @test norm(ru, Inf) <= 1e-12

    return nothing
end

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
θxs = copy(θ)
θ[uind] .+= θwz
θ .\= 1.0

F = cholesky(Symmetric(A*Diagonal(θ)*A'))

ξp = rand(m)
ξd = rand(n)
ξu = rand(1)

test_augmented_system(A, F, θ, θxs, θwz, uind, ξp, ξd, ξu)