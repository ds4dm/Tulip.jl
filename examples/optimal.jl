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
            0 <= x1, x2 <= 1
    =#
    m = TLP.Model{Tv}()
    m.params.OutputLevel = 1

    # Read problem and solve
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_opt.mps"))
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_Optimal

    # TODO: check validity of optimal solution
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_optimal(Float64)
end