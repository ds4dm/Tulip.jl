print("\tmodel.jl")

m = 8
n = 8
A = sprand(m, n, 0.5)
c = 1.0 - 2.0*rand(n)
b = 1.0 + rand(m)
u = 10.0 + sprand(n, 0.8)
p = nnz(u)

# Tests for PrimalDualPoint
v1 = Tulip.PrimalDualPoint(rand(n), rand(p), rand(m), rand(n), rand(p))
v2 = copy(v1)
@test n == sum(v1.x .== v2.x)
@test p == sum(v1.w .== v2.w)
@test m == sum(v1.y .== v2.y)
@test n == sum(v1.s .== v2.s)
@test p == sum(v1.z .== v2.z)

v3 = v1 + v2
@test n == sum(v3.x .== (v1.x + v2.x))
@test p == sum(v3.w .== (v1.w + v2.w))
@test m == sum(v3.y .== (v1.y + v2.y))
@test n == sum(v3.s .== (v1.s + v2.s))
@test p == sum(v3.z .== (v1.z + v2.z))


# Create random instance and run checks
model = Tulip.Model(A, b, c, u.nzind, u.nzval)

@test m == model.n_con
@test n == model.n_var
@test model.status == :Built

sol = model.sol
@test typeof(sol) <: Tulip.PrimalDualPoint
@test n == size(sol.x, 1)
@test p == size(sol.w, 1)
@test m == size(sol.y, 1)
@test n == size(sol.s, 1)
@test p == size(sol.z, 1)

println("\tPassed.")