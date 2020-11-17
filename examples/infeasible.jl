using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_infeasible(::Type{Tv};
    atol::Tv = 100 * sqrt(eps(Tv)),
    rtol::Tv = 100 * sqrt(eps(Tv)),
    kwargs...
) where{Tv}
    #=
    Infeasible example

    min     x1 + x2
    s.t.    x1 + x2 =  1
            x1 - x2 =  0
                 x2 =  1
            x1, x2  >= 0
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1
    # Set optional parameters
    for (k, val) in kwargs
        TLP.set_parameter(m, String(k), val)
    end

    # Read problem from .mps file and solve    
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_inf.mps"))
    TLP.optimize!(m)

    # Check status
    @test TLP.get_attribute(m, TLP.Status()) == TLP.Trm_PrimalInfeasible
    z = TLP.get_attribute(m, TLP.ObjectiveValue())
    @test z == Inf
    @test m.solution.primal_status == TLP.Sln_Unknown
    @test m.solution.dual_status   == TLP.Sln_InfeasibilityCertificate

    # Check unbounded dual ray
    y1 = m.solution.y_lower[1] - m.solution.y_upper[1]
    y2 = m.solution.y_lower[2] - m.solution.y_upper[2]
    y3 = m.solution.y_lower[3] - m.solution.y_upper[3]
    s1 = m.solution.s_lower[1] - m.solution.s_upper[1]
    s2 = m.solution.s_lower[2] - m.solution.s_upper[2]

    @test y1 + y3 >= atol  # dual cost should be > 0
    @test isapprox(y1 + y2      + s1, 0, atol=atol, rtol=rtol)
    @test isapprox(y1 - y2 + y3 + s2, 0, atol=atol, rtol=rtol)
    @test s1 >= -atol
    @test s2 >= -atol
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_infeasible(Float64)
end
