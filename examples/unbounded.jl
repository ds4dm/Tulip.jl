using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

function ex_unbounded(::Tv) where{Tv<:Real}
    #=
    Unbounded example

    min     -x1 + -x2
            x1 - x2 =  1
            x1, x2  >= 0
    =#
    m = TLP.Model{Tv}()

    # Create example
    x1 = TLP.add_variable!(m, "x1", -1.0, 0.0, Inf)
    x2 = TLP.add_variable!(m, "x2", -1.0, 0.0, Inf)
    c1 = TLP.add_constraint!(m, "c1", 1.0, 1.0, [x1, x2], [1,  -1])

    m.env.verbose = 1
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_DualInfeasible

    # TODO: check validity of infeasibility certificate
end

ex_unbounded(0.0)