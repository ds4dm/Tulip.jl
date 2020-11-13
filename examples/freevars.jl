using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_freevars(::Type{Tv};
    atol::Tv = 100 * sqrt(eps(Tv)),
    rtol::Tv = 100 * sqrt(eps(Tv)),
    kwargs...
) where{Tv}
    #=
    Example with free variables

    min       x1 +   x2 + x3
    s.t.    2 x1 +   x2      >= 2
              x1 + 2 x2      >= 2
              x1 +   x2 + x3 >= 0
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1
    # Set optional parameters
    for (k, val) in kwargs
        setfield!(m.params, k, val)
    end

    # Read problem and solve
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_freevars.mps"))
    TLP.optimize!(m)

    # Check status
    @test TLP.get_attribute(m, TLP.Status()) == TLP.Trm_Optimal
    z = TLP.get_attribute(m, TLP.ObjectiveValue())
    @test isapprox(z, 0, atol=atol, rtol=rtol)
    @test m.solution.primal_status == TLP.Sln_Optimal
    @test m.solution.dual_status   == TLP.Sln_Optimal

    # Check optimal solution
    x1 = m.solution.x[1]
    x2 = m.solution.x[2]
    x3 = m.solution.x[3]

    # Check primal feasibility (note there's no unique solution)
    @test 2*x1 +   x2      >= 2 - atol
    @test   x1 + 2*x2      >= 2 - atol
    @test   x1 +   x2 + x3 >= -atol
    
    # Free variables should have zero reduced cost
    s1 = m.solution.s_lower[1] - m.solution.s_upper[1]
    s2 = m.solution.s_lower[2] - m.solution.s_upper[2]
    s3 = m.solution.s_lower[3] - m.solution.s_upper[3]
    @test isapprox(s1, 0, atol=atol, rtol=rtol)
    @test isapprox(s2, 0, atol=atol, rtol=rtol)
    @test isapprox(s3, 0, atol=atol, rtol=rtol)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_freevars(Float64)
end