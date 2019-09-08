const MOI = Tulip.MOI
const MOIT = MOI.Test

const OPTIMIZER = TLP.Optimizer()

const CONFIG = MOIT.TestConfig()

@testset "Unit Tests" begin
    MOIT.basic_constraint_tests(OPTIMIZER, CONFIG; delete=false, get_constraint_function=false, get_constraint_set=false)
    # MOIT.unittest(OPTIMIZER, MOIT.TestConfig(atol=1e-6))
    # MOIT.modificationtest(OPTIMIZER, CONFIG)
end

# @testset "Linear tests" begin
#     @testset "Default Solver"  begin
#         MOIT.contlineartest(OPTIMIZER, MOIT.TestConfig(basis = true), [
#             # This requires an infeasiblity certificate for a variable bound.
#             "linear12"
#         ])
#     end
#     @testset "No certificate" begin
#         MOIT.linear12test(OPTIMIZER, MOIT.TestConfig(infeas_certificates=false))
#     end
# end