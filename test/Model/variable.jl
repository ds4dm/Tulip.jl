# TODO: run this with multiple combinations of Tv types.
# This may require to encapsulate each testset in its own typed function.
# Tv types: Float32, Float64, BigFloat, Rational
function run_tests_variable(::Tv) where{Tv<:Real}
    # Constructors
    @testset "VarId" begin
        vid = TLP.VarId(1)

        @test vid.uuid == 1
    end

    @testset "Bounds" begin

        # These should pass
        @test TLP._check_bounds(TLP.TLP_BND_UP, Tv(-Inf), Tv(0))
        @test TLP._check_bounds(TLP.TLP_BND_LO, Tv(0), Tv(Inf))
        @test TLP._check_bounds(TLP.TLP_BND_FX, Tv(1), Tv(1))
        @test TLP._check_bounds(TLP.TLP_BND_FR, Tv(-Inf), Tv(Inf))
        @test TLP._check_bounds(TLP.TLP_BND_RG, Tv(0), Tv(1))

        # These should fail
        # TODO: add all possible cases

    end

    @testset "VarData" begin
        vdat = TLP.VarData{Tv}("x1", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))

        @test vdat.name == "x1"
        @test vdat.obj == zero(Tv)
        @test vdat.bt == TLP.TLP_BND_LO
        @test vdat.lb == zero(Tv)
        @test vdat.ub == Tv(Inf)

        vdat = TLP.VarData("x1", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))

        @test isa(vdat, TLP.VarData{Tv})
        @test vdat.name == "x1"
        @test vdat.obj == zero(Tv)
        @test vdat.bt == TLP.TLP_BND_LO
        @test vdat.lb == zero(Tv)
        @test vdat.ub == Tv(Inf)
    end
        
    @testset "Constructors" begin
        vid = TLP.VarId(1)
        vdat = TLP.VarData{Tv}("x1", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))

        v = TLP.Variable(vid, vdat)
        @test isa(v, TLP.Variable{Tv})

        v = TLP.Variable{Tv}(vid, "x1", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
        @test isa(v, TLP.Variable{Tv})
    end

    vid = TLP.VarId(1)
    vdat = TLP.VarData("x1", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
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

        TLP.set_bounds!(v, TLP.TLP_BND_RG, Tv(-1), Tv(1))

        @test TLP.get_bounds(v) == (TLP.TLP_BND_RG, Tv(-1), Tv(1))

        @test TLP.get_lower_bound(v) == Tv(-1)
        @inferred TLP.get_lower_bound(v)
        @test TLP.get_upper_bound(v) ==  Tv(1)
        @inferred TLP.get_upper_bound(v)

        # Check type conversion
        TLP.set_bounds!(v, TLP.TLP_BND_RG, Tv(0), Tv(2))
        @test TLP.get_lower_bound(v) == Tv(0)
        @inferred TLP.get_lower_bound(v)

        @test TLP.get_upper_bound(v) ==  Tv(2)
        @inferred TLP.get_upper_bound(v)
    end

    return nothing
end

@testset "Variable" begin
    for Tv in [Float32, Float64, BigFloat, Rational]
        @testset "$Tv" begin run_tests_variable(zero(Tv)) end
    end
end 