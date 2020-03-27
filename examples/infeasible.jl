using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_infeasible(Tv::Type)
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

    # Read problem from .mps file and solve    
    TLP.load_problem!(m, joinpath(INSTANCE_DIR, "lpex_inf.mps"))
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_PrimalInfeasible

    # TODO: check validity of infeasibility certificate
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_infeasible(Float64)
end
