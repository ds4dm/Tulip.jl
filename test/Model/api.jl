# TODO: unit tests with various numerical types: Float32, Float64, BigFloat
function run_tests_api(::Tv) where{Tv<:Real}

    m = Tulip.Model_{Float64}()
    env = m.env

    x1 = TLP.add_variable!(m, "x1", 1.0, TLP.TLP_BND_LO, 0.0, Inf)
    x2 = TLP.add_variable!(m, "x1", 2.0, TLP.TLP_BND_LO, 0.0, Inf)

    c1 = TLP.add_constraint!(m, "c1", TLP.TLP_BND_FX, 2.0, 2.0, [x1, x2], [1.0, 1.0])

    # TODO: build standard form and solve model

    return nothing
end

@testset "low-level API" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_api(zero(Tv)) end
    end
end