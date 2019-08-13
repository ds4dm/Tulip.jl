import Tulip: HSDSolver, optimize!

using Random

function run_tests_hsd(::Tv) where{Tv<:Real}
    @testset "step length" begin
        m, n, p = 2, 2, 1
        pt = TLP.Point{Tv}(
            m, n, p,
            ones(Tv, n), ones(Tv, p), one(Tv),
            zeros(Tv, m), ones(Tv, n), ones(Tv, p), one(Tv),
            one(Tv)
        )
        d = TLP.Point{Tv}(
            m, n, p,
            -ones(Tv, n), -ones(Tv, p), -one(Tv),
            zeros(Tv, m), -ones(Tv, n), -ones(Tv, p), -one(Tv),
            one(Tv)
        )

        # Max step length for a single (x, d)
        @inferred TLP.max_step_length(pt.x, d.x)
        @test TLP.max_step_length(ones(Tv, 1), ones(Tv, 1)) == Tv(Inf)
        @test TLP.max_step_length(ones(Tv, 1), -ones(Tv, 1)) ≈ one(Tv)
        @test TLP.max_step_length(zeros(Tv, 1), -ones(Tv, 1)) ≈ zero(Tv)
        @test TLP.max_step_length(zeros(Tv, 1), ones(Tv, 1)) == Tv(Inf)

        # Max step length for the whole primal-dual point
        @inferred TLP.max_step_length(pt, d)
        @test TLP.max_step_length(pt, d) ≈ one(Tv)
    end

    @testset "Residuals" begin

    end

    @testset "Convergence" begin

        # TODO: optimal

        # TODO: dual infeasible

        # TODO: primal infeasible

        # TODO: ill-posed

    end
end

@testset "HSD" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_hsd(zero(Tv)) end
    end
end