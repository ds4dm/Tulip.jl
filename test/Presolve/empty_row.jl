function empty_row_tests(Tv::Type)

    # Build the following model
    #=
        min     x + y
        s.t.   -1 ⩽ 0 * x + 0 * y + 0 * z ⩽ 1
                1 ⩽ 0 * x + 0 * y + 0 * z ⩽ 2
    =#
    pb = Tulip.ProblemData{Tv}()

    m, n = 2, 3
    A = spzeros(Tv, m, n)

    b = ones(Tv, m)
    c = ones(Tv, n)

    Tulip.load_problem!(pb, "test",
        true, c, zero(Tv),
        A, Tv.([-1, 1]), Tv.([1, 2]), zeros(Tv, n), fill(Tv(Inf), n),
        ["c1", "c2"], ["x", "y", "z"]
    )

    ps = Tulip.PresolveData(pb)

    @test !ps.updated
    @test ps.nzrow[1] == ps.nzrow[2] == 0

    # Remove first empty row
    Tulip.remove_empty_row!(ps, 1)

    @test ps.updated
    @test ps.status == Tulip.Trm_Unknown
    @test ps.nrow == 1
    @test !ps.rowflag[1] && ps.rowflag[2]
    @test length(ps.ops) == 1

    op = ps.ops[1]
    @test isa(op, Tulip.EmptyRow{Tv})
    @test op.i == 1
    @test iszero(op.y)

    # Remove second empty row
    # This should detect infeasibility
    Tulip.remove_empty_row!(ps, 2)

    @test ps.status == Tulip.Trm_PrimalInfeasible
    @test ps.nrow == 1
    @test !ps.rowflag[1] && ps.rowflag[2]
    @test length(ps.ops) == 1

    # Check solution status & objective value
    sol = ps.solution
    @test sol.dual_status == Tulip.Sln_InfeasibilityCertificate
    @test sol.z_primal == sol.z_dual == Tv(Inf)

    # Check Farkas ray
    #   (current problem only has 1 row)
    @test sol.y_lower[1] >  zero(Tv)

    return
end

@testset "Empty row" begin
    for Tv in TvTYPES
        @testset "$Tv" begin empty_row_tests(Tv) end
    end
end