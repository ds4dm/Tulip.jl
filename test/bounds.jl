function run_tests_bounds(::Tv) where{Tv<:Real}

    # Valid bounds
    @test TLP.bound_type(Tv(-Inf), Tv(0)) == TLP.TLP_UP
    @test TLP.bound_type(Tv(0), Tv(Inf)) == TLP.TLP_LO
    @test TLP.bound_type(Tv(1), Tv(1)) == TLP.TLP_FX
    @test TLP.bound_type(Tv(-Inf), Tv(Inf)) == TLP.TLP_FR
    @test TLP.bound_type(Tv(0), Tv(1)) == TLP.TLP_RG

    # Invalid bounds
    @test_throws ErrorException TLP.bound_type(Tv(Inf), Tv(Inf))
    @test_throws ErrorException TLP.bound_type(Tv(1), Tv(0))
end

@testset "Bounds" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_bounds(zero(Tv)) end
    end
end