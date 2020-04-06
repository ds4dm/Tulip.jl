using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

INSTANCE_DIR = joinpath(@__DIR__, "dat")

function ex_unbounded(Tv::Type)
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
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_DualInfeasible

    # TODO: check validity of infeasibility certificate
end

if abspath(PROGRAM_FILE) == @__FILE__
    ex_unbounded(Float64)
end