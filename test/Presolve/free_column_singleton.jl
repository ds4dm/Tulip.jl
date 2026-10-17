function test_free_column_singleton_zero_dual(T::Type)
    # Free column singleton with zero cost and a negative coefficient (issue #163)
    #=
        min     x1
        s.t.    x1 - x2 ⩽ 0
                1 ⩽ x1 ⩽ 2,  x2 free
    =#
    # The dual of the row is y = 0 / -1 = -0.0 and the row lower bound is -Inf,
    # so selecting the bound from the sign of y gave -0.0 * -Inf = NaN.
    c, A = T.([1, 0]), sparse(T.([1 -1]))
    lr, ur = T.([-Inf]), T.([0])
    lc, uc = T.([1, -Inf]), T.([2, Inf])

    pb = Tulip.ProblemData{T}()
    Tulip.load_problem!(pb, "test", true, c, zero(T), A, lr, ur, lc, uc, ["r"], ["x1", "x2"])

    ps = Tulip.PresolveData(pb)
    Tulip.remove_free_column_singleton!(ps, 2)

    @test length(ps.ops) == 1
    @test isa(ps.ops[1], Tulip.FreeColumnSingleton{T})
    @test !ps.rowflag[1] && !ps.colflag[2]
    @test ps.obj0 == zero(T)

    # The objective value and the solution must be finite
    m = Tulip.Model{T}()
    Tulip.load_problem!(m.pbdata, "test", true, c, zero(T), A, lr, ur, lc, uc, ["r"], ["x1", "x2"])
    Tulip.set_parameter(m, "OutputLevel", 0)
    Tulip.optimize!(m)

    atol = 100 * sqrt(eps(T))
    @test Tulip.get_attribute(m, Tulip.Status()) == Tulip.Trm_Optimal
    @test isapprox(m.solution.z_primal, one(T); atol = atol)
    @test all(isfinite, m.solution.x)
    @test m.solution.x[1] - m.solution.x[2] <= atol

    return nothing
end

@testset "Free column singleton" begin
    for T in TvTYPES
        @testset "$T" begin
            test_free_column_singleton_zero_dual(T)
        end
    end
end
