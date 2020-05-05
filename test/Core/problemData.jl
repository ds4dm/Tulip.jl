# using Test

# include(joinpath(@__DIR__, "../../src/Core/problemData.jl"))

# TvTYPES = [Float64, BigFloat]

function run_tests_pbdata(Tv::Type)

    @testset "Creation" begin
        pb = TLP.ProblemData{Tv}("test")

        @test pb.name == "test"

        check_problem_size(pb, 0, 0)
        @test pb.objsense
        @test iszero(pb.obj0)

        # Add two columns
        #=
            min     x1 + 2 x2
            s.t.    0 ⩽ x1 ⩽ ∞
                    1 ⩽ x2 ⩽ ∞
        =#
        TLP.add_variable!(pb, Int[], Tv[],     one(Tv), zero(Tv), Tv(Inf), "x1")
        TLP.add_variable!(pb, Int[], Tv[], 2 * one(Tv),  one(Tv), Tv(Inf), "x2")
        
        check_problem_size(pb, 0, 2)
        col1, col2 = pb.acols[1], pb.acols[2]
        @test pb.obj == [one(Tv), 2*one(Tv)]
        @test pb.lvar == [zero(Tv), one(Tv)]
        @test pb.uvar == [Tv(Inf), Tv(Inf)]
        @test length(col1.nzind) == length(col1.nzval) == 0
        @test length(col2.nzind) == length(col2.nzval) == 0
        @test pb.var_names == ["x1", "x2"]

        # Add two constraints
        #=
            min     x1 + 2 x2
            s.t.    -∞ ⩽  -x1 +   x2 ⩽ 1
                    -1 ⩽ 2 x1 - 2 x2 ⩽ 0
                    0 ⩽ x1 ⩽ ∞
                    1 ⩽ x2 ⩽ ∞
        =#
        TLP.add_constraint!(pb, [1, 2], Tv.([-1, 1]), Tv(-Inf), one(Tv), "row1")
        TLP.add_constraint!(pb, [1, 2], Tv.([2, -2]), -one(Tv), zero(Tv), "row2")
        
        # Check dimensions
        check_problem_size(pb, 2, 2)

        # Check coefficients
        row1, row2 = pb.arows[1], pb.arows[2]
        @test row1.nzind == [1, 2]
        @test row1.nzval == Tv.([-1, 1])
        @test row2.nzind == [1, 2]
        @test row2.nzval == Tv.([2, -2])
        @test col1.nzind == [1, 2]
        @test col1.nzval == Tv.([-1, 2])
        @test col2.nzind == [1, 2]
        @test col2.nzval == Tv.([1, -2])

        # Check row bounds
        @test pb.lcon == [Tv(-Inf), -one(Tv)]
        @test pb.ucon == [one(Tv), zero(Tv)]
        # Check names
        @test pb.con_names == ["row1", "row2"]
        @test pb.var_names == ["x1", "x2"]


        # TODO: move these to HSD tests
        @testset "HSD" begin
            params = TLP.Parameters{Tv}()
            hsd = TLP.HSDSolver{Tv}(params, pb)

            A = hsd.kkt.A
            @test Matrix(A) == Tv.([
                [-1 1 1 0];
                [2 -2 0 -1]
            ])

            params.OutputLevel = 0
            TLP.optimize!(hsd, params)

            @test hsd.primal_bound_scaled ≈ Tv(5 // 2)
        end  # HSD test

        empty!(pb)
        @test pb.name == ""
        @test iszero(pb.obj0)
        check_problem_size(pb, 0, 0)

    end

    @testset "Modification" begin

    end

    @testset "Queries" begin

    end

    return nothing
end

function check_problem_size(pb::TLP.ProblemData, ncon::Int, nvar::Int)
    @test pb.ncon == ncon
    @test pb.nvar == nvar

    @test length(pb.obj) == nvar

    @test length(pb.arows) == ncon
    @test length(pb.acols) == nvar

    @test length(pb.lcon) == ncon
    @test length(pb.ucon) == ncon
    @test length(pb.lvar) == nvar
    @test length(pb.uvar) == nvar

    @test length(pb.con_names) == ncon
    @test length(pb.var_names) == nvar
    return nothing
end

@testset "ProblemData" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_pbdata(Tv) end
    end
end