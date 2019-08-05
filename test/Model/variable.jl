using Test

include("../../src/Model/variable.jl")

# TODO: run this with multiple combinations of (Tv, Ti) types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
# Ti types: Int32, Int64, Int128, BigInt
function run_tests_variable()
    # Constructors
    @testset "Constructors" begin
        vid = VarId(1)
        vdat = VarData("x1", 0.0, 0.0, Inf)
        v = Variable(vid, vdat)

        v = Variable{Float64, Int}(vid)

        v = Variable{Float64, Int}(vid, "x1", 0.0, 0.0, Inf)
    end

    vid = VarId(1)
    vdat = VarData("x1", 0.0, 0.0, Inf)
    v = Variable(vid, vdat)

    # UUID
    @testset "UUID" begin
        @test get_uuid(v) == vid
    end

    # Name
    @testset "Name" begin
        @test get_name(v) == "x1"

        set_name!(v, "x2")
        @test get_name(v) == "x2"
    end

    # Objective coefficients
    @testset "Objective" begin
        @test get_obj_coeff(v) == 0.0

        # With initial type
        set_obj_coeff!(v, 1.0)
        @test get_obj_coeff(v) == 1.0
        @inferred get_obj_coeff(v)

        # Check type conversion
        set_obj_coeff!(v, 0)
        @inferred get_obj_coeff(v)
        @test get_obj_coeff(v) == 0.0
    end

    # Bounds
    @testset "Bounds" begin
        @test get_lower_bound(v) == 0.0
        @test get_upper_bound(v) == Inf

        set_lower_bound!(v, -1.0)
        set_upper_bound!(v,  1.0)

        @test get_lower_bound(v) == -1.0
        @inferred get_lower_bound(v)
        @test get_upper_bound(v) ==  1.0
        @inferred get_upper_bound(v)

        # Check type conversion
        set_lower_bound!(v, 0)
        @test get_lower_bound(v) == 0.0
        @inferred get_lower_bound(v)

        set_upper_bound!(v, 2)
        @test get_upper_bound(v) ==  2.0
        @inferred get_upper_bound(v)
    end

    return nothing
end

@testset "Variable" begin run_tests_variable() end 