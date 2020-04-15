using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_unbounded(::Type{Tv};
    atol::Tv = 100 * sqrt(eps(Tv)),
    rtol::Tv = 100 * sqrt(eps(Tv))
) where{Tv}
    #=
    Unbounded example

    min     -x1 + -x2
            x1 - x2 =  1
            x1, x2  >= 0
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1

    # Read problem from .mps file and solve    
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_ubd.mps"))
    TLP.optimize!(m)

    # Check status
    @test TLP.get_attribute(m, TLP.Status()) == TLP.Trm_DualInfeasible
    z = TLP.get_attribute(m, TLP.ObjectiveValue())
    @test z == -Tv(Inf)
    @test m.solution.primal_status == TLP.Sln_InfeasibilityCertificate
    @test m.solution.dual_status   == TLP.Sln_Unknown

    # Check unbounded ray
    x1 = m.solution.x[1]
    x2 = m.solution.x[2]
    Ax1 = m.solution.Ax[1]

    @test x1 >= -atol
    @test x2 >= -atol
    @test isapprox(Ax1, 0, atol=atol, rtol=rtol)
    @test -x1 - x2 <= -atol
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_unbounded(Float64)
end