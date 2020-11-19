using Test, MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities

const OPTIMIZER = TLP.Optimizer()

MOI.set(OPTIMIZER, MOI.Silent(), true)

const CONFIG = MOIT.TestConfig(basis=false, atol=1e-6, rtol=1e-6)

const MOI_EXCLUDE = [
    # Unit tests
    "delete_nonnegative_variables",         # Requires Vector-Of-Variables
    "update_dimension_nonnegative_variables",
    "delete_soc_variables",
    "solve_func_vectoraffine_nonneg",       # Requires `VectorAffineFunction` and `Nonnegatives`
    # Require SingleVariable objective function
    "get_objective_function",
    "solve_result_index",
    "solve_singlevariable_obj",
    "solve_single_variable_dual_max",
    "solve_single_variable_dual_min",
    "solve_integer_edge_cases",             # Requires integer variables
    "solve_qcp_edge_cases",                 # Requires quadratic constraints
    "solve_qp_edge_cases",                  # Requires quadratic objective
    "solve_zero_one_with_bounds_1",         # Requires binary variables
    "solve_zero_one_with_bounds_2",         # Requires binary variables
    "solve_zero_one_with_bounds_3",         # Requires binary variables
    "solve_affine_deletion_edge_cases",     # Requires VectorAffineFunction-in-Nonpositives
    "solve_objbound_edge_cases",            # Requires integer variables
    "solve_duplicate_terms_vector_affine",  # Requires `VectorAffineFunction`
    "variablenames",                        # TODO, requires to settle the convention on variable names
    "raw_status_string",                    # TODO
    "solve_time",                           # TODO
    # Modifications
    "delete_variable_with_single_variable_obj",  # Requires SingleVariable objective function
    "solve_const_vectoraffine_nonpos",      # Requires VectorAffineFunction-in-Nonpositives
    "solve_multirow_vectoraffine_nonpos",   # Requires VectorAffineFunction-in-Nonpositives
    # Linear tests
    "linear7",                              # Requires `VectorAffineFunction` and `Nonnegatives`/`Nonpositives`
    "linear15",                             # Requires `VectorAffineFunction` and `Zeros`
    "partial_start"                         # Requires `VariablePrimalStart`
]

@testset "MOI Unit Tests" begin
    @testset "Basic constraint tests" begin
        MOIT.basic_constraint_tests(OPTIMIZER, CONFIG)
    end

    @testset "Other unit tests" begin
        MOIT.unittest(OPTIMIZER, CONFIG, MOI_EXCLUDE)
    end

    @testset "Modifications" begin
        MOIT.modificationtest(OPTIMIZER, CONFIG, MOI_EXCLUDE)
    end
end

# Run the MOI tests with HSD and MPC algorithms
for ipm in [Tulip.HSD, Tulip.MPC]
    @testset "MOI Linear tests - $ipm" begin
        OPTIMIZER.inner.params.IPM.Factory = Tulip.Factory(ipm)
        MOIT.contlineartest(OPTIMIZER, CONFIG, MOI_EXCLUDE)
    end
end

MOI.empty!(OPTIMIZER)
const BRIDGED = MOI.Bridges.full_bridge_optimizer(OPTIMIZER, Float64)
MOI.set(BRIDGED, MOI.Silent(), true)

@testset "Bridged tests" begin
    MOIT.unittest(BRIDGED, CONFIG, [
        "get_objective_function",               # Requires SingleVariable objective
        "delete_soc_variables",                 # Requires SOC constraints
        "solve_integer_edge_cases",             # Requires integer variables
        "solve_qcp_edge_cases",                 # Requires quadratic constraints
        "solve_qp_edge_cases",                  # Requires quadratic objective
        "solve_zero_one_with_bounds_1",         # Requires binary variables
        "solve_zero_one_with_bounds_2",         # Requires binary variables
        "solve_zero_one_with_bounds_3",         # Requires binary variables
        "solve_objbound_edge_cases",            # Requires integer variables
        "variablenames",                        # TODO, requires to settle the convention on variable names
        "raw_status_string",                    # TODO
        "solve_time",                           # TODO
    ])
    MOIT.contlineartest(BRIDGED, CONFIG, ["partial_start"])
end
