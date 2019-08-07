# TODO: tests in multiple precision: Float32, Float64, BigFloat
function run_tests_standardform()

    model = TLP.Model_{Float64}()

    A0 = [
        [1 -1 1 1   -1   1 0  0  0.0];
        [0  0 0 4.4 -4.4 0 1  0  0.0];
        [0 -2 0 0    0   0 0 -1  0.0];
        [0  0 0 0    0   5 0  0 -1.0]
    ]

    # Add each kind of variable
    x_fx = TLP.add_variable!(model, "xfx", 1.0, TLP.TLP_BND_FX,  0.1, 0.1)
    x_up = TLP.add_variable!(model, "xup", 2.0, TLP.TLP_BND_UP, -Inf, 0.02)
    x_lo = TLP.add_variable!(model, "xlo", 3.0, TLP.TLP_BND_LO,  0.003, Inf)
    x_fr = TLP.add_variable!(model, "xfr", 4.0, TLP.TLP_BND_FR, -Inf, Inf)
    x_rg = TLP.add_variable!(model, "xrg", 5.0, TLP.TLP_BND_RG,  0.000050, 0.000055)

    # Add each kind of constraint
    c_fx = TLP.add_constraint!(model, "cfx", TLP.TLP_BND_FX, 1.0, 1.0,
        # This row will go from
        # [1, 1, 1, 1, 1] to
        # [1, -1, 1, 1, -1, 1, 0, 0, 0]
        [x_fx, x_up, x_lo, x_fr, x_rg], ones(5)
    )
    c_up = TLP.add_constraint!(model, "cup", TLP.TLP_BND_UP, -Inf, 20.0,
        # This row will go from
        # [0, 0, 0, 4,  0] to
        # [0, 0, 0, 4, -4, 0, 1, 0, 0]
        [x_fr], [4.4]
    )
    c_lo = TLP.add_constraint!(model, "clo", TLP.TLP_BND_LO, 300.0, Inf,
        # This row will go from
        # [0,  2, 0, 0, 0] to
        # [0, -2, 0, 0, 0, 0, 0, -1, 0]
        [x_up], [2.0]
    )
    c_rg = TLP.add_constraint!(model, "crg", TLP.TLP_BND_RG, -4000.0, 4000.0,
        # This row will go from
        # [0, 0, 0, 0, 5] to
        # [0, 0, 0, 0, 0, 5, 0, 0, -1]
        [x_rg], [5.0]
    )
    # free constraints are not allowed

    # Convert to standard form
    ncons, nvars, aI, aJ, aV, b, c, uind, uval, con2idx, var2idx, idx2con, idx2var = TLP.convert_to_standard_form(model.pbdata_raw)

    @test ncons == 4
    @test nvars == 9  # 5 vars + 1 free split + 3 slacks

    # objective
    @test c[1] == 1.0
    @test c[2] == -2.0
    @test c[3] == 3.0
    @test c[4] == 4.0
    @test c[5] == -4.0
    @test c[6] == 5.0

    # Right-hand-side
    b0 = [1.0, 20.0, 300.0, -4000.0]
    l0 = [0.1, 0.02, 0.003, 0.0, 0.0, 0.000050, 0.0, 0.0, 0.0]

    @test b ≈ (b0 .- A0 * l0)

    # Constraint matrix
    @test length(aI) == length(aJ) == length(aV) == 13
    A = sparse(aI, aJ, aV, ncons, nvars)
    @test Matrix{Float64}(A) ≈ A0

    # Upper-bounds
    @test uind == [1, 6, 9]
    @test uval ≈ [0.0, 5e-6, 8000.0]

end

@testset "StandardForm" begin run_tests_standardform() end

# Simple example:
#=
    min     -x1 -x2
    s.t.    x1 + x2 in [0, 1]


model = TLP.Model_{Float64}()

x = TLP.add_variable!(model, "x", -1.0, TLP.TLP_BND_LO, 0.0, Inf)
y = TLP.add_variable!(model, "y", -1.0, TLP.TLP_BND_LO, 0.0, Inf)

c = TLP.add_constraint!(model, "c1", TLP.TLP_BND_RG, 0.0, 1.0, [x, y], [1.0, 1.0])

m, n, aI, aJ, aV, b, c, uind, uval, con2idx, var2idx, idx2con, idx2var = TLP.convert_to_standard_form(model.pbdata_raw)

@show m, n

A = sparse(aI, aJ, aV, m, n)
@show Matrix(A)

@show b
@show c

@show uind
@show uval
=#