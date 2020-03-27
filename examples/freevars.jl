using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_freevars(Tv::Type)
    #=
    Example with free variables

    min       x1 +   x2 + x3
    s.t.    2 x1 +   x2      >= 2
              x1 + 2 x2      >= 2
              x1 +   x2 + x3 >= 0
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1

    # Read problem and solve
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_freevars.mps"))
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_Optimal
    @test isapprox(hsd.primal_bound_scaled, 0, atol = sqrt(eps(Tv)))

    # TODO: Check optimal solution
end

ex_freevars(Float64)