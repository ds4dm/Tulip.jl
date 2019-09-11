const MOI = Tulip.MOI
const MOIT = MOI.Test

const OPTIMIZER = TLP.Optimizer()

const CONFIG = MOIT.TestConfig()

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
    ])
    # MOIT.modificationtest(OPTIMIZER, CONFIG)
end

# @testset "Linear tests" begin
#     @testset "Default Solver"  begin
#         MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis = false), String[
#             # This requires an infeasiblity certificate for a variable bound.
#             # "linear12"
#         ])
#     end
#     @testset "No certificate" begin
#         MOIT.linear12test(OPTIMIZER, MOIT.TestConfig(infeas_certificates=false))
#     end
# end