using Test

include("../../src/Model/constraint.jl")

# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
function run_tests_lincon()
    # Constructors
    @testset "Constructors" begin
        cid = ConstrId(1)
        cdat = LinConstrData{Float64}("c1", 0.0, Inf)
        c = LinearConstraint(cid, cdat)

        c = LinearConstraint{Float64}(cid)

        c = LinearConstraint{Float64}(cid, "c1", 0.0, Inf)
    end

    cid = ConstrId(1)
    cdat = LinConstrData{Float64}("c1", 0.0, Inf)
    c = LinearConstraint(cid, cdat)

    # UUID
    @testset "UUID" begin
        @test get_uuid(c) == cid
    end

    # Name
    @testset "Name" begin
        @test get_name(c) == "c1"

        set_name!(c, "c2")
        @test get_name(c) == "c2"
    end

    # Bounds
    @testset "Bounds" begin
        @test get_lower_bound(c) == 0.0
        @test get_upper_bound(c) == Inf

        set_lower_bound!(c, -1.0)
        set_upper_bound!(c,  1.0)

        @test get_lower_bound(c) == -1.0
        @inferred get_lower_bound(c)
        @test get_upper_bound(c) ==  1.0
        @inferred get_upper_bound(c)

        # Check type conversion
        set_lower_bound!(c, 0)
        @test get_lower_bound(c) == 0.0
        @inferred get_lower_bound(c)

        set_upper_bound!(c, 2)
        @test get_upper_bound(c) ==  2.0
        @inferred get_upper_bound(c)
    end

    return nothing
end

@testset "LinearConstraint" begin run_tests_lincon() end 