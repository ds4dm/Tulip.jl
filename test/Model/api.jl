# TODO: unit tests with various numerical types: Float32, Float64, BigFloat
function run_tests_api(::Tv) where{Tv<:Real}

    m = Tulip.Model_{Tv}()
    env = m.env

    @testset "add_variable" begin
        m = TLP.Model_{Tv}()

        @test TLP.get_num_var(m) == 0

        x = TLP.add_variable!(m)
        @test TLP.get_num_var(m) == 1

        y = TLP.add_variable!(m, "y", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
        @test TLP.get_num_var(m) == 2

        z = TLP.add_variable!(m, "z", 0, TLP.TLP_BND_FX, 1, 1, TLP.ConstrId[], Int[])
        @test TLP.get_num_var(m) == 3
        
        t = TLP.add_variable!(m, "t", zero(Tv), TLP.TLP_BND_RG, zero(Tv), oneunit(Tv), TLP.ConstrId[], Tv[])
        @test TLP.get_num_var(m) == 4
    end

    @testset "add_variables" begin
        # TODO
        #=
        m = TLP.Model_{Tv}()

        @test TLP.get_num_var(m) == 0
        x = TLP.add_variables!(m, 2)

        @test TLP.get_num_var(m) == 2
        =#
    end

    @testset "names" begin
        m = TLP.Model_{Tv}()

        x = TLP.add_variable!(m)
        @test TLP.get_var_name(m, x) == ""

        y = TLP.add_variable!(m, "y", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
        @test TLP.get_var_name(m, y) == "y"

        z = TLP.add_variable!(m, "z", 0, TLP.TLP_BND_FX, 1, 1, TLP.ConstrId[], Int[])
        @test TLP.get_var_name(m, z) == "z"
        
        t = TLP.add_variable!(m, "t", zero(Tv), TLP.TLP_BND_RG, zero(Tv), oneunit(Tv), TLP.ConstrId[], Tv[])
        @test TLP.get_var_name(m, t) == "t"
    end

    @testset "Bounds" begin
        m = TLP.Model_{Tv}()

        # Variables are added ot the model with Float64 input to check for type conversion
        xup = TLP.add_variable!(m, "xup", zero(Tv), TLP.TLP_BND_UP, -Inf, 0.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xup)
        @test bt == TLP.TLP_BND_UP
        @test lb == Tv(-Inf)
        @test ub == zero(Tv)

        xlo = TLP.add_variable!(m, "xlo", zero(Tv), TLP.TLP_BND_LO, 0.0, Inf)
        (bt, lb, ub) = TLP.get_var_bounds(m, xlo)
        @test bt == TLP.TLP_BND_LO
        @test lb == zero(Tv)
        @test ub == Tv(Inf)

        xfx = TLP.add_variable!(m, "xfx", zero(Tv), TLP.TLP_BND_FX, 1.0, 1.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xfx)
        @test bt == TLP.TLP_BND_FX
        @test lb == Tv(1.0)
        @test ub == Tv(1.0)

        xfr = TLP.add_variable!(m, "xfr", zero(Tv), TLP.TLP_BND_FR, -Inf, Inf)
        (bt, lb, ub) = TLP.get_var_bounds(m, xfr)
        @test bt == TLP.TLP_BND_FR
        @test lb == Tv(-Inf)
        @test ub == Tv(Inf)

        xrg = TLP.add_variable!(m, "xrg", zero(Tv), TLP.TLP_BND_RG, 0.0, 1.0)
        (bt, lb, ub) = TLP.get_var_bounds(m, xrg)
        @test bt == TLP.TLP_BND_RG
        @test lb == zero(Tv)
        @test ub == Tv(1.0)

    end

    x1 = TLP.add_variable!(m, "x1", 1.0, TLP.TLP_BND_LO, 0.0, Inf)
    x2 = TLP.add_variable!(m, "x1", 2.0, TLP.TLP_BND_RG, 0.0, 1.0)

    c1 = TLP.add_constraint!(m, "c1", TLP.TLP_BND_FX, 2.0, 2.0, [x1, x2], [1.0, 1.0])

    # TODO: build standard form and solve model
    # Do no include tests here yet
    # These will be functional tests for later
    #=
    ncons, nvars, aI, aJ, aV, b, c, uind, uval, con2idx, var2idx, idx2con, idx2var = TLP.convert_to_standard_form(m.pbdata_raw)

    A = sparse(aI, aJ, aV, ncons, nvars)
    F = TLP.symbolic_cholesky(A)

    m.pbdata_std = TLP.StandardForm(A, b, c, uind, uval)

    pt = TLP.Solvers.Point(
        ncons, nvars, length(uind),
        zeros(Tv, nvars), zeros(Tv, length(uind)), zero(Tv),
        zeros(Tv, ncons), zeros(Tv, nvars), zeros(Tv, length(uind)), zero(Tv), zero(Tv)
    )
    res = TLP.Solvers.Residuals(
        zeros(Tv, ncons), zeros(Tv, length(uind)), zeros(Tv, nvars), zero(Tv),
        zero(Tv), zero(Tv), zero(Tv), zero(Tv)
    )

    hsd = TLP.Solvers.HSDSolver(m.pbdata_std,
        0, TLP.TerminationStatus(0), TLP.SolutionStatus(0), TLP.SolutionStatus(0),
        pt, res, F
    )

    env[:verbose] = 1
    env[:barrier_iter_max] = 10

    @time TLP.Solvers.optimize!(hsd, env)
    =#

    return nothing
end

# Sparse factorizations are not supported for eltypes other than Float64
# Until then, two solutions:
#   1. Convert to dense and use dense factorizations
#   2. Drop support for eltypes other than Float64
@testset "low-level API" begin
    for Tv in [Float64]
        @testset "$Tv" begin run_tests_api(zero(Tv)) end
    end
end