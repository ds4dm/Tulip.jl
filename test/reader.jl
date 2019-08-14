INSTANCE_DIR = joinpath(@__DIR__, "../dat/dummy/")
function run_tests_mps_reader(::Tv) where{Tv<:Real}
    m = TLP.Model{Tv}()

    TLP.readmps!(m, joinpath(INSTANCE_DIR, "lpex_opt.mps"))

    @testset "Name" begin
        @test m.name == "LP1"
    end

    @testset "Variables" begin
        @test TLP.get_num_var(m) == 2

        vars = collect(keys(m.pbdata_raw.vars))
        @test length(vars) == 2
        
        x1, x2 = vars[1], vars[2]

        # Names
        @testset "Names" begin
            @test TLP.get_var_name(m, x1) == "X1"
            @test TLP.get_var_name(m, x2) == "X2"
        end

        # Objective coeffs
        @testset "Obj" begin
            @test m.pbdata_raw.vars[x1].dat.obj == Tv(1.0)
            @test m.pbdata_raw.vars[x2].dat.obj == Tv(2.0)
        end

        # Bounds
        @testset "Bounds" begin
            @test TLP.get_var_bounds(m, x1) == (TLP.TLP_RG, Tv(0.0), Tv(1.0))
            @test TLP.get_var_bounds(m, x2) == (TLP.TLP_RG, Tv(0.0), Tv(1.0))
        end
    end

    @testset "Constraints" begin
        @test TLP.get_num_constr(m.pbdata_raw) == 2

        cons = collect(keys(m.pbdata_raw.constrs))
        @test length(cons) ==2

        c1, c2 = cons[1], cons[2]

        @testset "Names" begin
            @test TLP.get_name(m.pbdata_raw.constrs[c1]) == "ROW1"
            @test TLP.get_name(m.pbdata_raw.constrs[c2]) == "ROW2"
        end

        @testset "Bounds" begin
            @test TLP.get_bounds(m.pbdata_raw.constrs[c1]) == (TLP.TLP_FX, Tv(1.0), Tv(1.0))
            @test TLP.get_bounds(m.pbdata_raw.constrs[c2]) == (TLP.TLP_FX, Tv(0.0), Tv(0.0))
        end
    end
end

# TODO: per-function unit tests
@testset "MPS reader" begin 
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_mps_reader(zero(Tv)) end
    end
end