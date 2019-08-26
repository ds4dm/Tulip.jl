function run_tests_mps_reader(::Tv) where{Tv<:Real}
    m = TLP.Model{Tv}()

    TLP.readmps!(m, joinpath(@__DIR__, "test_lp.mps"))

    @testset "Name" begin
        @test m.name == "TESTLP"
    end

    @testset "Variables" begin
        @test TLP.get_num_var(m) == 9

        x = collect(keys(m.pbdata_raw.vars))
        @test length(x) == 9

        # Names
        @testset "Names" begin
            for j in 1:9
                @test TLP.get_var_name(m, x[j]) == "X$(j)"
            end
        end

        # Objective coeffs
        @testset "Obj" begin
            for j in 1:9
                @test m.pbdata_raw.vars[x[j]].dat.obj == Tv(j)
            end
        end

        # Bounds
        @testset "Bounds" begin
            @test TLP.get_var_bounds(m, x[1]) == (TLP.TLP_LO, Tv(1), Tv(Inf))
            @test TLP.get_var_bounds(m, x[2]) == (TLP.TLP_RG, Tv(0), Tv(1))
            @test TLP.get_var_bounds(m, x[3]) == (TLP.TLP_FX, Tv(1), Tv(1))
            @test TLP.get_var_bounds(m, x[4]) == (TLP.TLP_FR, Tv(-Inf), Tv(Inf))
            @test TLP.get_var_bounds(m, x[5]) == (TLP.TLP_UP, Tv(-Inf), Tv(5))
            @test TLP.get_var_bounds(m, x[6]) == (TLP.TLP_LO, Tv(0), Tv(Inf))
            @test TLP.get_var_bounds(m, x[7]) == (TLP.TLP_RG, Tv(0), Tv(1))
            @test TLP.get_var_bounds(m, x[8]) == (TLP.TLP_LO, Tv(1), Tv(Inf))
            @test TLP.get_var_bounds(m, x[9]) == (TLP.TLP_RG, Tv(0), Tv(1))
        end
    end

    @testset "Constraints" begin
        @test TLP.get_num_constr(m.pbdata_raw) == 3

        cons = collect(keys(m.pbdata_raw.constrs))
        @test length(cons) == 3


        @testset "Names" begin
            for i in 1:3
                @test TLP.get_name(m.pbdata_raw.constrs[cons[i]]) == "ROW$(i)"
            end
        end

        @testset "Bounds" begin
            @test TLP.get_bounds(m.pbdata_raw.constrs[cons[1]]) == (TLP.TLP_FX, Tv(1), Tv(1))
            @test TLP.get_bounds(m.pbdata_raw.constrs[cons[2]]) == (TLP.TLP_UP, Tv(-Inf), Tv(2))
            @test TLP.get_bounds(m.pbdata_raw.constrs[cons[3]]) == (TLP.TLP_LO, Tv(3), Tv(Inf))
        end
    end
end

# TODO: per-function unit tests
@testset "MPS reader" begin 
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_mps_reader(zero(Tv)) end
    end
end