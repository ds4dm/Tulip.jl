using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_optimal(Tv::Type)
    #=
    Bounded example

    min     x1 + 2*x2
    s.t.    x1 + x2 =  1
            x1 - x2 =  0
            0 <= x1 <= 1
            0 <= x2 <= 1
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1

    # Read problem and solve
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_opt.mps"))
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_Optimal
    @test TLP.get_attribute(m, TLP.ObjectiveValue()) ≈ 1.5

    # Check validity of optimal solution
    x1 = m.solution.x[1]
    x2 = m.solution.x[2]
    Ax1 = m.solution.Ax[1]
    Ax2 = m.solution.Ax[2]

    @test isapprox(x1, 1 // 2, rtol = 100 * sqrt(eps(Tv)))
    @test isapprox(x2, 1 // 2, rtol = 100 * sqrt(eps(Tv)))
    @test isapprox(Ax1, 1)
    @test isapprox(Ax2, 0, atol = 100 * sqrt(eps(Tv)))

    y1 = m.solution.y_lower[1] - m.solution.y_upper[1]
    y2 = m.solution.y_lower[2] - m.solution.y_upper[2]
    s1 = m.solution.s_lower[1] - m.solution.s_upper[1]
    s2 = m.solution.s_lower[2] - m.solution.s_upper[2]

    @test isapprox(y1,  3 // 2, atol = 100 * sqrt(eps(Tv)))
    @test isapprox(y2, -1 // 2, atol = 100 * sqrt(eps(Tv)))
    @test isapprox(s1, 0, atol = 100 * sqrt(eps(Tv)))
    @test isapprox(s2, 0, atol = 100 * sqrt(eps(Tv)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_optimal(Float64)
end