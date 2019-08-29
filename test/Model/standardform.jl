# TODO: tests in multiple precision: Float32, Float64, BigFloat
function run_tests_standardform()

    model = TLP.Model{Float64}()

    A0 = [
        [1 -1 1 1   -1   1 0  0  0.0];
        [0  0 0 4.4 -4.4 0 1  0  0.0];
        [0 -2 0 0    0   0 0 -1  0.0];
        [0  0 0 0    0   5 0  0 -1.0]
    ]

    # Add each kind of variable
    x_fx = TLP.add_variable!(model, "xfx", 1.0,      0.1,      0.1)
    x_up = TLP.add_variable!(model, "xup", 2.0,     -Inf,     0.02)
    x_lo = TLP.add_variable!(model, "xlo", 3.0,    0.003,      Inf)
    x_fr = TLP.add_variable!(model, "xfr", 4.0,     -Inf,      Inf)
    x_rg = TLP.add_variable!(model, "xrg", 5.0, 0.000050, 0.000055)

    # Add each kind of constraint
    c_fx = TLP.add_constraint!(model, "cfx", 1.0, 1.0,
        # This row will go from
        # [1, 1, 1, 1, 1] to
        # [1, -1, 1, 1, -1, 1, 0, 0, 0]
        [x_fx, x_up, x_lo, x_fr, x_rg], ones(5)
    )
    c_up = TLP.add_constraint!(model, "cup", -Inf, 20.0,
        # This row will go from
        # [0, 0, 0, 4,  0] to
        # [0, 0, 0, 4, -4, 0, 1, 0, 0]
        [x_fr], [4.4]
    )
    c_lo = TLP.add_constraint!(model, "clo", 300.0, Inf,
        # This row will go from
        # [0,  2, 0, 0, 0] to
        # [0, -2, 0, 0, 0, 0, 0, -1, 0]
        [x_up], [2.0]
    )
    c_rg = TLP.add_constraint!(model, "crg", -4000.0, 4000.0,
        # This row will go from
        # [0, 0, 0, 0, 5] to
        # [0, 0, 0, 0, 0, 5, 0, 0, -1]
        [x_rg], [5.0]
    )
    # free constraints are not allowed

    # Convert to standard form
    sf = TLP.convert_to_standard_form(SparseMatrixCSC, model.pbdata_raw)

    @test sf.ncon == 4
    @test sf.nvar == 9  # 5 vars + 1 free split + 3 slacks
    @test sf.nupb == 3  # 1 FX var + 1 RG var + 1 RG con

    # objective
    @test sf.c[1] == 1.0
    @test sf.c[2] == -2.0
    @test sf.c[3] == 3.0
    @test sf.c[4] == 4.0
    @test sf.c[5] == -4.0
    @test sf.c[6] == 5.0

    # Right-hand-side
    b0 = [1.0, 20.0, 300.0, -4000.0]
    l0 = [0.1, 0.02, 0.003, 0.0, 0.0, 0.000050, 0.0, 0.0, 0.0]

    @test sf.b ≈ (b0 .- A0 * l0)

    # Constraint matrix
    @test Matrix{Float64}(sf.A) ≈ A0

    # Upper-bounds
    @test sf.uind == [1, 6, 9]
    @test sf.uval ≈ [0.0, 5e-6, 8000.0]

end

@testset "StandardForm" begin run_tests_standardform() end