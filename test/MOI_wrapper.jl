const MOI = Tulip.MOI
const MOIT = MOI.Test

const OPTIMIZER = TLP.Optimizer()

const CONFIG = MOIT.TestConfig(basis=false)

@testset "Unit Tests" begin
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG; delete=false)
    MOIT.unittest(OPTIMIZER, MOIT.TestConfig(atol=1e-6), [
        "solve_integer_edge_cases",         # Requires integer variables
        "solve_qcp_edge_cases",             # Requires quadratic constraints
        "solve_qp_edge_cases",              # Requires quadratic objective
        "solve_zero_one_with_bounds_1",     # Requires binary variables
        "solve_zero_one_with_bounds_2",     # Requires binary variables
        "solve_zero_one_with_bounds_3",     # Requires binary variables
        "solve_affine_deletion_edge_cases", # Requires VectorAffineFunction-in-Nonpositives
        "solve_objbound_edge_cases",         # Requires integer variables
        "solve_duplicate_terms_vector_affine",  # Requires `VectorAffineFunction`
    ])
    MOIT.modificationtest(OPTIMIZER, CONFIG)
end

@testset "Linear tests" begin
    MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis=false, atol=1e-6, rtol=1e-6), String[
        "linear7"   # Requires `VectorAffineFunction` and `Nonnegatives`/`Nonpositives`
        "linear15"  # Requires `VectorAffineFunction` and `Zeros`
        "partial_start_test"  # Requires `VariablePrimalStart`
    ])
end