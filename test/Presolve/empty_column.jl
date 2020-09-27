function emtpy_column_tests(T::Type)

    # We test all the following combinations:
    #=
        min          c * x
        s.t.    lb ≤     x ≤ ub


        ------------------------------
                    |        c       
        (lb, ub)    | +1  | -1  |  0
        ------------------------------
        (-∞,  u)    | -∞  |  u  |  u 
        (-∞, +∞)    | -∞  | +∞  |  0 
        ( l,  u)    |  l  |  u  |  l 
        ( l, +∞)    |  l  | +∞  |  l 
        ------------------------------
    =#
    function build_problem(l, u, c)
        pb = Tulip.ProblemData{T}()

        Tulip.load_problem!(pb, "Test", 
            true, [c], zero(T),
            spzeros(T, 0, 1), T[], T[], [l], [u], String[], ["x"]
        )
        return pb
    end

    L = T.([-Inf,  -1])
    U = T.([   1, Inf])
    C = T.([-1, 0, 1])
    for l in L, u in U, c in C
        @testset "$((l, u, c))" begin
            pb = build_problem(l, u, c)

            ps = Tulip.PresolveData(pb)

            # Remove empty variable
            Tulip.remove_empty_column!(ps, 1)

            if c > 0 && !isfinite(l)
                @test ps.status == Tulip.Trm_DualInfeasible
                @test ps.colflag[1]
                @test ps.ncol == 1

                sol = ps.solution
                @test sol.primal_status == Tulip.Sln_InfeasibilityCertificate
                @test sol.m == 0 && sol.n == 1
                @test sol.x[1] < 0
            elseif c < 0 && !isfinite(u)
                @test ps.status == Tulip.Trm_DualInfeasible
                @test ps.colflag[1]
                @test ps.ncol == 1

                sol = ps.solution
                @test sol.primal_status == Tulip.Sln_InfeasibilityCertificate
                @test sol.m == 0 && sol.n == 1
                @test sol.x[1] > 0
            else
                @test ps.status == Tulip.Trm_Unknown
                @test !ps.colflag[1]
                @test ps.ncol == 0
                @test ps.updated

                # Check that operation was recorded correctly
                @test length(ps.ops) == 1
                op = ps.ops[1]
                @test isa(op, Tulip.EmptyColumn)
                @test op.j == 1
            end 
        end  # testset
    end  # loop

    return
end

@testset "Empty column" begin
    for T in TvTYPES
        @testset "$T" begin emtpy_column_tests(T) end
    end
end