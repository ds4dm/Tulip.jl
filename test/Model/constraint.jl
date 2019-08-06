# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
function run_tests_lincon()
    # Constructors
    @testset "Constructors" begin
        cid = TLP.ConstrId(1)
        cdat = TLP.LinConstrData{Float64}("c1", TLP.TLP_BND_LO, 0.0, Inf)
        c = TLP.LinearConstraint(cid, cdat)

        c = TLP.LinearConstraint{Float64}(cid)

        c = TLP.LinearConstraint{Float64}(cid, "c1", TLP.TLP_BND_LO, 0.0, Inf)
    end

    cid = TLP.ConstrId(1)
    cdat = TLP.LinConstrData{Float64}("c1", TLP.TLP_BND_LO, 0.0, Inf)
    c = TLP.LinearConstraint(cid, cdat)

    # UUID
    @testset "UUID" begin
        @test TLP.get_uuid(c) == cid
    end

    # Name
    @testset "Name" begin
        @test TLP.get_name(c) == "c1"

        TLP.set_name!(c, "c2")
        @test TLP.get_name(c) == "c2"
    end

    # Bounds
    @testset "Bounds" begin
        @test TLP.get_lower_bound(c) == 0.0
        @test TLP.get_upper_bound(c) == Inf

        TLP.set_bounds!(c, TLP.TLP_BND_RG, -1.0, 1.0)

        @test TLP.get_lower_bound(c) == -1.0
        @inferred TLP.get_lower_bound(c)
        @test TLP.get_upper_bound(c) ==  1.0
        @inferred TLP.get_upper_bound(c)

        # Check type conversion
        TLP.set_bounds!(c, TLP.TLP_BND_RG, 0, 2)
        @test TLP.get_lower_bound(c) == 0.0
        @inferred TLP.get_lower_bound(c)

        @test TLP.get_upper_bound(c) ==  2.0
        @inferred TLP.get_upper_bound(c)
    end

    return nothing
end

@testset "LinearConstraint" begin run_tests_lincon() end 