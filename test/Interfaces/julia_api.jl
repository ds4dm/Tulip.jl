import Tulip
using Test

function test_reader()
    lp = Tulip.Model{Float64}()

    Tulip.load_problem!(lp, joinpath(@__DIR__, "lp.mps"))
    check_data(lp)

    Tulip.load_problem!(lp, joinpath(@__DIR__, "lp.mps.gz"))
    check_data(lp)

    Tulip.load_problem!(lp, joinpath(@__DIR__, "lp.mps.bz2"))
    check_data(lp)

end

function check_data(lp)

    pb = lp.pbdata

    @test pb.name == "LP1"

    @test pb.ncon == 2
    @test pb.nvar == 2

    @test pb.objsense
    @test pb.obj0 == 0.0
    @test pb.obj == [1.0, 2.0]

    @test pb.lvar == [0.0, 0.0]
    @test pb.uvar == [1.0, 1.0]
    @test pb.lcon == [1.0, 0.0]
    @test pb.ucon == [1.0, 0.0]

    @test pb.con_names == ["ROW1", "ROW2"]
    @test pb.var_names == ["X1", "X2"]

    col1, col2 = pb.acols[1], pb.acols[2]
    @test col1.nzind == [1, 2]
    @test col1.nzval == [1.0, 1.0]
    @test col2.nzind == [1, 2]
    @test col2.nzval == [1.0, -1.0]

    row1, row2 = pb.arows[1], pb.arows[2]
    @test row1.nzind == [1, 2]
    @test row1.nzval == [1.0, 1.0]
    @test row2.nzind == [1, 2]
    @test row2.nzval == [1.0, -1.0]



end

@testset "Reader" begin
    test_reader()
end
