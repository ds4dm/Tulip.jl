# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat
function run_tests_lincon(::Tv) where{Tv<:Real}
    
    @testset "ConstrId" begin
        cid = TLP.ConstrId(1)

        @test cid.uuid == 1
    end

    @testset "LinConstrData" begin
        cdat = TLP.LinConstrData{Tv}("c1", zero(Tv), oneunit(Tv))

        @test cdat.name == "c1"
        @test cdat.bt == TLP.TLP_RG
        @test cdat.lb == zero(Tv)
        @test cdat.ub == oneunit(Tv)
    end

    @testset "Constructors" begin
        cid = TLP.ConstrId(1)
        cdat = TLP.LinConstrData{Tv}("c1", zero(Tv), Tv(Inf))

        # Construct from ID and (typed) data
        c = TLP.LinearConstraint(cid, cdat)
        @test isa(c, TLP.LinearConstraint{Tv})

        # Construct from ID and (typed) data fields
        c = TLP.LinearConstraint{Tv}(cid, "c1", zero(Tv), Tv(Inf))
        @test isa(c, TLP.LinearConstraint{Tv})

        # Construct from ID and un-typed data fields
        c = TLP.LinearConstraint{Tv}(cid, "c1", 0, Inf)
        @test isa(c, TLP.LinearConstraint{Tv})
    end

    cid = TLP.ConstrId(1)
    cdat = TLP.LinConstrData{Tv}("c1", zero(Tv), Tv(Inf))
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

        TLP.set_bounds!(c, Tv(-1), Tv(1))

        @inferred TLP.get_bounds(c)
        @test TLP.get_bounds(c) == (TLP.TLP_RG, Tv(-1), Tv(1))

        # Individual getters for bounds
        @inferred TLP.get_lower_bound(c)
        @inferred TLP.get_upper_bound(c)
        @test TLP.get_lower_bound(c) == Tv(-1)
        @test TLP.get_upper_bound(c) ==  Tv(1)
        
        # Check type conversion
        TLP.set_bounds!(c, 0, 2)
        
        @inferred TLP.get_lower_bound(c)
        @inferred TLP.get_upper_bound(c)
        @test TLP.get_upper_bound(c) ==  Tv(2)
        @test TLP.get_lower_bound(c) == zero(Tv)
    end

    return nothing
end

@testset "LinearConstraint" begin 
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_lincon(zero(Tv)) end
    end
end 