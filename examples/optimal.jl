using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "../dat/dummy")

function ex_optimal(::Tv) where{Tv<:Real}
    #=
    Bounded example

    min     x1 + 2*x2
    s.t.    x1 + x2 =  1
            x1 - x2 =  0
            0 <= x1, x2 <= 1
    =#
    m = TLP.Model{Tv}()
    m.env.verbose = 1

    # Read problem and solve
    TLP.readmps!(m, joinpath(INSTANCE_DIR, "lpex_opt.mps"))
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_Optimal

    # TODO: check validity of optimal solution
end

ex_optimal(0.0)