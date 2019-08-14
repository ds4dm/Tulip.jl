using LinearAlgebra
using SparseArrays
using Test

import Tulip
TLP = Tulip

function ex_optimal(::Tv) where{Tv<:Real}
    #=
    Bounded example

    min     x1 + x2
    s.t.    x1 + x2 =  1
            x1 - x2 =  0
            0 <= x1, x2 <= 1
    =#
    m = TLP.Model{Tv}()

    # Create example
    x1 = TLP.add_variable!(m, "x1", 1.0, 0.0, 1.0)
    x2 = TLP.add_variable!(m, "x2", 1.0, 0.0, 1.0)
    c1 = TLP.add_constraint!(m, "c1", 1.0, 1.0, [x1, x2], [1,  1])
    c2 = TLP.add_constraint!(m, "c2", 0.0, 0.0, [x1, x2], [1, -1])

    m.env.verbose = 1
    TLP.optimize!(m)

    # Check status
    hsd = m.solver
    @test hsd.solver_status == TLP.Trm_Optimal

    # TODO: check validity of optimal solution
end

ex_optimal(0.0)