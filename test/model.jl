print("\tmodel.jl")

m = 8
n = 8
A = sprand(m, n, 0.5)
c = 1.0 - 2.0*rand(n)
b = 1.0 + rand(m)
u = 10.0 + sprand(n, 0.8)
p = nnz(u)

# Create random instance and run checks
model = Tulip.Model(A, b, c, u.nzind, u.nzval)

@test m == model.n_con
@test n == model.n_var
@test model.status == :Built

@test n == size(model.x, 1)
@test p == size(model.w, 1)
@test m == size(model.y, 1)
@test n == size(model.s, 1)
@test p == size(model.z, 1)

println("\tPassed.")