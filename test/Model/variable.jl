# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat, Rational
function run_tests_variable(::Tv) where{Tv<:Real}
    
    @testset "VarId" begin
        vid = TLP.VarId(1)

        @test vid.uuid == 1
    end

    @testset "VarData" begin

        # internal constructor
        vdat = TLP.VarData{Tv}("x1", zero(Tv), zero(Tv), Tv(Inf))
        @test vdat.name == "x1"
        @test vdat.obj == zero(Tv)
        @test vdat.bt == TLP.TLP_LO
        @test vdat.lb == zero(Tv)
        @test vdat.ub == Tv(Inf)

        # Constructor with un-typed data fields
        vdat = TLP.VarData{Tv}("x1", 0, 0, Inf)
        @test isa(vdat, TLP.VarData{Tv})
        @test vdat.name == "x1"
        @test vdat.obj == zero(Tv)
        @test vdat.bt == TLP.TLP_LO
        @test vdat.lb == zero(Tv)
        @test vdat.ub == Tv(Inf)
    end
        
    @testset "Constructors" begin
        vid = TLP.VarId(1)
        vdat = TLP.VarData{Tv}("x1", zero(Tv), zero(Tv), Tv(Inf))

        # Construct from ID and (typed) data
        v = TLP.Variable(vid, vdat)
        @test isa(v, TLP.Variable{Tv})

        # Construct from ID and (typed) data fields
        v = TLP.Variable{Tv}(vid, "x1", zero(Tv), zero(Tv), Tv(Inf))
        @test isa(v, TLP.Variable{Tv})

        # Construct from ID an un-typed data fields
        v = TLP.Variable{Tv}(vid, "x1", 0, 0, Inf)
        @test isa(v, TLP.Variable{Tv})
    end

    vid = TLP.VarId(1)
    vdat = TLP.VarData{Tv}("x1", zero(Tv), zero(Tv), Tv(Inf))
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
        @test TLP.get_obj_coeff(v) == zero(Tv)

        # With initial type
        TLP.set_obj_coeff!(v, Tv(1))
        @test TLP.get_obj_coeff(v) == Tv(1)
        @inferred TLP.get_obj_coeff(v)

        # Check type conversion
        TLP.set_obj_coeff!(v, zero(Tv))
        @inferred TLP.get_obj_coeff(v)
        @test TLP.get_obj_coeff(v) == zero(Tv)
    end

    # Bounds
    @testset "Bounds" begin
        @test TLP.get_lower_bound(v) == zero(Tv)
        @test TLP.get_upper_bound(v) == Tv(Inf)

        TLP.set_bounds!(v, Tv(-1), Tv(1))

        @inferred TLP.get_bounds(v)
        @test TLP.get_bounds(v) == (TLP.TLP_RG, Tv(-1), Tv(1))

        # Check for type stability
        @inferred TLP.get_lower_bound(v)
        @inferred TLP.get_upper_bound(v)

        # check values
        @test TLP.get_lower_bound(v) == Tv(-1)
        @test TLP.get_upper_bound(v) == Tv(1)

        # Check type conversion
        TLP.set_bounds!(v, 0, 2)

        @inferred TLP.get_lower_bound(v)
        @inferred TLP.get_upper_bound(v)
        @test TLP.get_lower_bound(v) == Tv(0)
        @test TLP.get_upper_bound(v) == Tv(2)
        
    end

    return nothing
end

@testset "Variable" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_variable(zero(Tv)) end
    end
end 