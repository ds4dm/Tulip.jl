# TODO: unit tests with various numerical types: Float32, Float64, BigFloat
function run_tests_api(::Tv) where{Tv<:Real}

    m = Tulip.Model{Tv}()
    env = m.env

    @testset "add_variable" begin
        m = TLP.Model{Tv}()

        @test TLP.get_num_var(m) == 0

        x = TLP.add_variable!(m)
        @test TLP.get_num_var(m) == 1

        y = TLP.add_variable!(m, "y", zero(Tv), zero(Tv), Tv(Inf))
        @test TLP.get_num_var(m) == 2
        @test TLP.get_var_name(m, y) == "y"

        z = TLP.add_variable!(m, "z", 0, 1, 1, TLP.ConstrId[], Int[])
        @test TLP.get_num_var(m) == 3
        @test TLP.get_var_name(m, z) == "z"
        
        t = TLP.add_variable!(m, "t", zero(Tv), zero(Tv), oneunit(Tv), TLP.ConstrId[], Tv[])
        @test TLP.get_num_var(m) == 4
        @test TLP.get_var_name(m, t) == "t"
    end

    @testset "add_variables" begin
        # TODO
        #=
        m = TLP.Model{Tv}()

        @test TLP.get_num_var(m) == 0
        x = TLP.add_variables!(m, 2)

        @test TLP.get_num_var(m) == 2
        =#
    end

    @testset "Bounds" begin
        m = TLP.Model{Tv}()

        # Variables are added ot the model with Float64 input to check for type conversion
        xup = TLP.add_variable!(m, "xup", zero(Tv), -Inf, 0.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xup)
        @test bt == TLP.TLP_UP
        @test lb == Tv(-Inf)
        @test ub == zero(Tv)

        xlo = TLP.add_variable!(m, "xlo", zero(Tv), 0.0, Inf)
        (bt, lb, ub) = TLP.get_var_bounds(m, xlo)
        @test bt == TLP.TLP_LO
        @test lb == zero(Tv)
        @test ub == Tv(Inf)

        xfx = TLP.add_variable!(m, "xfx", zero(Tv), 1.0, 1.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xfx)
        @test bt == TLP.TLP_FX
        @test lb == Tv(1.0)
        @test ub == Tv(1.0)

        xfr = TLP.add_variable!(m, "xfr", zero(Tv), -Inf, Inf)
        (bt, lb, ub) = TLP.get_var_bounds(m, xfr)
        @test bt == TLP.TLP_FR
        @test lb == Tv(-Inf)
        @test ub == Tv(Inf)

        xrg = TLP.add_variable!(m, "xrg", zero(Tv), 0.0, 1.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xrg)
        @test bt == TLP.TLP_RG
        @test lb == zero(Tv)
        @test ub == Tv(1.0)

    end

    return nothing
end

# Sparse factorizations are not supported for eltypes other than Float64
# Until then, two solutions:
#   1. Convert to dense and use dense factorizations
#   2. Drop support for eltypes other than Float64
#   3. Use the LDLFactorization package of Dominique
@testset "low-level API" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_api(zero(Tv)) end
    end
end