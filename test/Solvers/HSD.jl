import Tulip.Solvers:
    HSDSolver, optimize!

using Random

function test_hsd()
    m, n = 4, 8

    Random.seed!(0)
    A = sprand(m, n, 1.0)
    b = A * ones(n)
    c = rand(n)
    uind = Int.(collect(1:n))
    uval = 10.0 .* ones(n)

    pb = TLP.StandardForm(A, b, c, uind, uval)
    env = Tulip.TulipEnv()

    # Instantiate HSD solver
    hsd = HSDSolver{Float64}(pb)

    # solve
    optimize!(hsd, env)

    # Check status
    # println(hsd.solver_status)
    # println(hsd.primal_status)
    # println(hsd.dual_status)
end

@testset "HSD" begin
    test_hsd()
end