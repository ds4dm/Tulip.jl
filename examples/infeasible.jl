using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

function ex_infeasible(::Tv) where{Tv<:Real}
    #=
    Infeasible example

    min     x1 + x2
    s.t.    x1 + x2 =  1
            x1 - x2 =  0
                 x2 =  1
            x1, x2  >= 0
    =#
    m = TLP.Model{Tv}()

    # Create example
    x1 = TLP.add_variable!(m, "x1", 1.0, 0.0, Inf)
    x2 = TLP.add_variable!(m, "x2", 1.0, 0.0, Inf)
    c1 = TLP.add_constraint!(m, "c1", 1.0, 1.0, [x1, x2], [1,  1])
    c2 = TLP.add_constraint!(m, "c2", 0.0, 0.0, [x1, x2], [1, -1])
    c3 = TLP.add_constraint!(m, "c3", 1.0, 1.0, [x1, x2], [0,  1])

    m.env.verbose = 1
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_PrimalInfeasible

    # TODO: check validity of infeasibility certificate
end

ex_infeasible(0.0)