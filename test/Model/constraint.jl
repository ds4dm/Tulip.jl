# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
function run_tests_lincon(::Tv) where{Tv<:Real}
    
    @testset "ConstrId" begin
        cid = TLP.ConstrId(1)

        @test cid.uuid == 1
    end

    @testset "LinConstrData" begin
        cdat = TLP.LinConstrData{Tv}("c1", TLP.TLP_BND_RG, zero(Tv), oneunit(Tv))

        @test cdat.name == "c1"
        @test cdat.bt == TLP.TLP_BND_RG
        @test cdat.lb == zero(Tv)
        @test cdat.ub == oneunit(Tv)
    end

    @testset "Constructors" begin
        cid = TLP.ConstrId(1)
        cdat = TLP.LinConstrData{Tv}("c1", TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
        c = TLP.LinearConstraint(cid, cdat)

        c = TLP.LinearConstraint{Tv}(cid, "c1", TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
    end

    cid = TLP.ConstrId(1)
    cdat = TLP.LinConstrData{Tv}("c1", TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
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
        @test TLP.get_lower_bound(c) == zero(Tv)
        @test TLP.get_upper_bound(c) == Tv(Inf)

        TLP.set_bounds!(c, TLP.TLP_BND_RG, Tv(-1), Tv(1))

        @test TLP.get_lower_bound(c) == Tv(-1)
        @inferred TLP.get_lower_bound(c)
        @test TLP.get_upper_bound(c) ==  Tv(1)
        @inferred TLP.get_upper_bound(c)

        # Check type conversion
        TLP.set_bounds!(c, TLP.TLP_BND_RG, 0, 2)
        @test TLP.get_lower_bound(c) == zero(Tv)
        @inferred TLP.get_lower_bound(c)

        @test TLP.get_upper_bound(c) ==  Tv(2)
        @inferred TLP.get_upper_bound(c)
    end

    return nothing
end

@testset "LinearConstraint" begin 
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_lincon(zero(Tv)) end
    end
end 