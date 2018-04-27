include("../src/ipm.jl")
include("../src/readmps.jl")

# read file
print("Reading file...")
m, n, obj, rhs, coeffs, lb, ub, ranges = readmps("../dat/netlib/AFIRO.SIF");
println(" Done.\n")

# solve model
tol = 10.0^-8
println("Solving model...")
x, y, s, X, Y, S = solve(coeffs, rhs, obj, tol, 1);

println()
println("Optimal objective value: ", @sprintf("%12.6e", dot(x, obj)))