function run_tests_env(::Tv) where{Tv<:Real}
    env = Tulip.Env{Tv}()

    env_ = copy(env)
    @test isa(env_, Tulip.Env{Tv})
    @test env_.verbose == env.verbose

end  # testset

@testset "Env" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_env(zero(Tv)) end
    end
end