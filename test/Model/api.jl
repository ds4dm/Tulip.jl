# TODO: unit tests with various numerical types: Float32, Float64, BigFloat
function run_tests_api()

    @testset "Model creation" begin
        # Create variable
        m = TLP.Model_("test", TLP.ProblemData{Float64}())

        # Create two variables
        x = TLP.add_variable!(m)
        @test x.uuid == 1

        y = TLP.add_variable!(m, "y", 1.0, TLP.TLP_BND_LO, 0.0, Inf)
        @test y.uuid == 2

        # Create a constraint
        c = TLP.add_constraint!(m, "c1", TLP.TLP_BND_FX, 1.0, 1.0, [x, y], [1.0, 2.0])
        @test c.uuid == 1

        # Check that coefficients were updated properly
        @test m.pbdata_raw.coeffs[x, c] == 1.0
        @test m.pbdata_raw.coeffs[y, c] == 2.0
    end

    return nothing
end

@testset "low-level API" begin run_tests_api() end