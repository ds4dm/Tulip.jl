using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_optimal(::Type{Tv};
    atol::Tv = 100 * sqrt(eps(Tv)),
    rtol::Tv = 100 * sqrt(eps(Tv)),
    kwargs...
) where{Tv}

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
    # Set optional parameters
    for (k, val) in kwargs
        setfield!(m.params, k, val)
    end

    # Read problem and solve
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_opt.mps"))
    TLP.optimize!(m)

    # Check status
    @test TLP.get_attribute(m, TLP.Status()) == TLP.Trm_Optimal
    z = TLP.get_attribute(m, TLP.ObjectiveValue())
    @test isapprox(z, 3 // 2, atol=atol, rtol=rtol)
    @test m.solution.primal_status == TLP.Sln_Optimal
    @test m.solution.dual_status   == TLP.Sln_Optimal

    # Check primal solution
    x1 = m.solution.x[1]
    x2 = m.solution.x[2]
    Ax1 = m.solution.Ax[1]
    Ax2 = m.solution.Ax[2]

    @test isapprox(x1, 1 // 2, atol=atol, rtol=rtol)
    @test isapprox(x2, 1 // 2, atol=atol, rtol=rtol)
    @test isapprox(Ax1, 1, atol=atol, rtol=rtol)
    @test isapprox(Ax2, 0, atol=atol, rtol=rtol)

    # Check duals
    y1 = m.solution.y_lower[1] - m.solution.y_upper[1]
    y2 = m.solution.y_lower[2] - m.solution.y_upper[2]
    s1 = m.solution.s_lower[1] - m.solution.s_upper[1]
    s2 = m.solution.s_lower[2] - m.solution.s_upper[2]
    @test isapprox(y1,  3 // 2, atol=atol, rtol=rtol)
    @test isapprox(y2, -1 // 2, atol=atol, rtol=rtol)
    @test isapprox(s1, 0, atol=atol, rtol=rtol)
    @test isapprox(s2, 0, atol=atol, rtol=rtol)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_optimal(Float64)
end