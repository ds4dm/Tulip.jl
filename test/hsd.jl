function test_augmented_system(A, F, θ, θxs, θwz, uind, ξp, ξd, ξu)

    dx = zeros(length(ξd))
    dy = zeros(length(ξp))
    dz = zeros(length(ξu))

    # Solve augmented system
    Tulip.solve_augsys_hsd!(A, F, θ, θwz, uind, dx, dy, dz, ξp, ξd, ξu)

    # Compute residuals
    rp = A*dx .- ξp

    dz_ = zeros(length(ξd))
    dz_[uind] .= dz
    rd = - θxs .* dx .+ transpose(A)*dy .- dz_ .- ξd

    ru = dx[uind] .- θwz .\ dz - ξu

    # Check numerics
    @test norm(rp, Inf) <= 1e-12
    @test norm(rd, Inf) <= 1e-12
    @test norm(ru, Inf) <= 1e-12

    return nothing
end

function test_newton_system(
    A, F, b, c, uind, uval, θ, θwz,
    p, q, r, ρ,
    x, w, y, s, z, t, k,
    ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk
)

    dx = zeros(length(ξd))
    dw = zeros(length(ξu))
    dy = zeros(length(ξp))
    ds = zeros(length(ξd))
    dz = zeros(length(ξu))
    dt = Ref(0.0)
    dk = Ref(0.0)

    # Solve augmented system
    Tulip.solve_newton_hsd!(
        A, F, b, c, uind, uval, θ, θwz,
        p, q, r, ρ,
        x, w, y, s, z, t, k,
        dx, dw, dy, ds, dz, dt, dk,
        ξp, ξu, ξd, ξg, ξxs, ξwz, ξtk
    )

    return nothing
end

function test_step_length()
    @test Tulip.max_step_length([1.0], [1.0]) == Inf
    @test Tulip.max_step_length([1.0], [-2.0]) == 0.5

    # This should raise a DimensionMismatch
    @test_throws DimensionMismatch Tulip.max_step_length([1.0, 1.0], [1.0])
    return nothing
end

m, n = 2, 3
Random.seed!(0)
A = sparse([
    [1.0 2.0 0.9];
    [0.4 0.5 0.6]
])
uind = [2]
uval = [4.0]

x = [0.2, 0.3, 0.4]
w = [0.5]
y = [0.9, -3.0]
s = [1.0, 2.0, 0.1]
z = [0.3]
t = Ref(1.0)
k = Ref(1.0)


θwz = w ./ z
θ = s ./ x
θxs = copy(θ)
θ[uind] .+= θwz
θ .\= 1.0

F = cholesky(Symmetric(A*Diagonal(θ)*A'))

ξp = [2.0, -1.0]
ξd = [1.0, 1.1, 0.9]
ξu = [0.1]
ξg = 0.1
ξxs = [0.9, 0.8, 0.5]
ξwz = [1.1]
ξtk = 4.0

test_augmented_system(A, F, θ, θxs, θwz, uind, ξp, ξd, ξu)

p = zeros(n)
q = zeros(m)
r = zeros(1)
b = [1.1, 0.9]
c = [9.9, 8.7, -0.4]

Tulip.solve_augsys_hsd!(
    A, F, θ, θwz, uind,
    p, q, r,
    b, c, uval
)
ρ = (k.x / t.x) - dot(c, p) + dot(b, q) - dot(uval, r)

test_newton_system(
    A, F, b, c, uind, uval, θ, θwz,
    p, q, r, ρ,
    x, w, y, s, z, t, k,
    ξp, ξd, ξu, ξg, ξxs, ξwz, ξtk
)