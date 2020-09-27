function run_tests_pbdata(T::Type)

    @testset "Creation" begin
        pb = TLP.ProblemData{T}("test")

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
        TLP.add_variable!(pb, Int[], T[],     one(T), zero(T), T(Inf), "x1")
        TLP.add_variable!(pb, Int[], T[], 2 * one(T),  one(T), T(Inf), "x2")
        
        check_problem_size(pb, 0, 2)
        col1, col2 = pb.acols[1], pb.acols[2]
        @test pb.obj == [one(T), 2*one(T)]
        @test pb.lvar == [zero(T), one(T)]
        @test pb.uvar == [T(Inf), T(Inf)]
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
        TLP.add_constraint!(pb, [1, 2], T.([-1, 1]), T(-Inf), one(T), "row1")
        TLP.add_constraint!(pb, [1, 2], T.([2, -2]), -one(T), zero(T), "row2")
        
        # Check dimensions
        check_problem_size(pb, 2, 2)

        # Check coefficients
        row1, row2 = pb.arows[1], pb.arows[2]
        @test row1.nzind == [1, 2]
        @test row1.nzval == T.([-1, 1])
        @test row2.nzind == [1, 2]
        @test row2.nzval == T.([2, -2])
        @test col1.nzind == [1, 2]
        @test col1.nzval == T.([-1, 2])
        @test col2.nzind == [1, 2]
        @test col2.nzval == T.([1, -2])

        # Check row bounds
        @test pb.lcon == [T(-Inf), -one(T)]
        @test pb.ucon == [one(T), zero(T)]
        # Check names
        @test pb.con_names == ["row1", "row2"]
        @test pb.var_names == ["x1", "x2"]

        empty!(pb)
        @test pb.name == ""
        @test iszero(pb.obj0)
        check_problem_size(pb, 0, 0)
    end

    @testset "Delete" begin
        pb = TLP.ProblemData{T}("test")
        #=
            min     x1 + 2 x2 + 3 x3
            s.t.    1 ⩽ 1 * x1 ⩽ 10
                    2 ⩽ 2 * x2 ⩽ 20
                    3 ⩽ 3 * x3 ⩽ 30

                    11 ⩽ x1 ⩽ 110
                    22 ⩽ x2 ⩽ 220
                    33 ⩽ x3 ⩽ 330
        =#

        TLP.add_variable!(pb, Int[], T[],     one(T), 11 * one(T), 110 * one(T), "x1")
        TLP.add_variable!(pb, Int[], T[], 2 * one(T), 22 * one(T), 220 * one(T), "x2")
        TLP.add_variable!(pb, Int[], T[], 3 * one(T), 33 * one(T), 330 * one(T), "x3")

        TLP.add_constraint!(pb, [1], T.([1]), 1 * one(T), 10 * one(T), "row1")
        TLP.add_constraint!(pb, [2], T.([2]), 2 * one(T), 20 * one(T), "row2")
        TLP.add_constraint!(pb, [3], T.([3]), 3 * one(T), 30 * one(T), "row3")

        # Delete row 1 and check remaining problem
        TLP.delete_constraint!(pb, 1)

        @test pb.ncon == 2
        @test pb.nvar == 3
        row2, row3 = pb.arows
        @test pb.con_names == ["row2", "row3"]
        @test pb.lcon == T.([2, 3])
        @test pb.ucon == T.([20, 30])

        @test row2.nzind == [2]
        @test row2.nzval == [T(2)]

        @test row3.nzind == [3]
        @test row3.nzval == [T(3)]

        # Delete variable 2
        TLP.delete_variable!(pb, 2)
        @test pb.ncon == 2
        @test pb.nvar == 2
        col1, col3 = pb.acols
        @test pb.var_names == ["x1", "x3"]
        @test pb.lvar == T.([11, 33])
        @test pb.uvar == T.([110, 330])

        @test col1.nzind == []
        @test col1.nzval == T[]

        @test col3.nzind == [2]
        @test col3.nzval == [T(3)]
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
    for T in TvTYPES
        @testset "$T" begin run_tests_pbdata(T) end
    end
end