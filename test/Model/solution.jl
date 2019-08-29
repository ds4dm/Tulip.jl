"""
    run_tests_solution

"""
function run_tests_solution(::Tv) where{Tv<:Real}

    @testset "UP" begin
        m = TLP.Model{Tv}()
        x = TLP.add_variable!(m, "x", -1.0, -Inf, 1.0)
        TLP.optimize!(m)

        @test TLP.get_value(m, x) ≈ Tv(1)
    end

    @testset "LO" begin
        m = TLP.Model{Tv}()
        x = TLP.add_variable!(m, "x", 1.0, 1.0, Inf)
        TLP.optimize!(m)

        @test TLP.get_value(m, x) ≈ Tv(1)
    end

    @testset "FX" begin
        m = TLP.Model{Tv}()
        x = TLP.add_variable!(m, "x", 1, 1, 1)
        TLP.optimize!(m)

        @test TLP.get_value(m, x) ≈ Tv(1)
    end

    @testset "FR" begin
        m = TLP.Model{Tv}()
        x = TLP.add_variable!(m, "x", -1, -Inf, Inf)
        con = TLP.add_constraint!(m, "con", 0.0, 1.0, [x], [1.0])
        TLP.optimize!(m)

        @test TLP.get_value(m, x) ≈ Tv(1)
    end

    @testset "RG" begin
        m = TLP.Model{Tv}()
        x = TLP.add_variable!(m, "x", -1, 0.0, 1.0)
        TLP.optimize!(m)

        @test TLP.get_value(m, x) ≈ Tv(1)
    end


end

@testset "Solution API" begin
    for Tv in [Float64]
        @testset "$Tv" begin run_tests_solution(zero(Tv)) end
    end
end