# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
function run_tests_variable()
    # Constructors
    @testset "Constructors" begin
        vid = TLP.VarId(1)
        vdat = TLP.VarData("x1", 0.0, TLP.TLP_BND_LO, 0.0, Inf)
        v = TLP.Variable(vid, vdat)

        v = TLP.Variable{Float64}(vid)

        v = TLP.Variable{Float64}(vid, "x1", 0.0, TLP.TLP_BND_LO, 0.0, Inf)
    end

    vid = TLP.VarId(1)
    vdat = TLP.VarData("x1", 0.0, TLP.TLP_BND_LO, 0.0, Inf)
    v = TLP.Variable(vid, vdat)

    # UUID
    @testset "UUID" begin
        @test TLP.get_uuid(v) == vid
    end

    # Name
    @testset "Name" begin
        @test TLP.get_name(v) == "x1"

        TLP.set_name!(v, "x2")
        @test TLP.get_name(v) == "x2"
    end

    # Objective coefficients
    @testset "Objective" begin
        @test TLP.get_obj_coeff(v) == 0.0

        # With initial type
        TLP.set_obj_coeff!(v, 1.0)
        @test TLP.get_obj_coeff(v) == 1.0
        @inferred TLP.get_obj_coeff(v)

        # Check type conversion
        TLP.set_obj_coeff!(v, 0)
        @inferred TLP.get_obj_coeff(v)
        @test TLP.get_obj_coeff(v) == 0.0
    end

    # Bounds
    @testset "Bounds" begin
        @test TLP.get_lower_bound(v) == 0.0
        @test TLP.get_upper_bound(v) == Inf

        TLP.set_bounds!(v, TLP.TLP_BND_RG, -1.0, 1.0)

        @test TLP.get_lower_bound(v) == -1.0
        @inferred TLP.get_lower_bound(v)
        @test TLP.get_upper_bound(v) ==  1.0
        @inferred TLP.get_upper_bound(v)

        # Check type conversion
        TLP.set_bounds!(v, TLP.TLP_BND_RG, 0, 2)
        @test TLP.get_lower_bound(v) == 0.0
        @inferred TLP.get_lower_bound(v)

        @test TLP.get_upper_bound(v) ==  2.0
        @inferred TLP.get_upper_bound(v)
    end

    return nothing
end

@testset "Variable" begin run_tests_variable() end 