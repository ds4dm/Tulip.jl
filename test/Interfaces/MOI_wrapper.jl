import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

const OPTIMIZER = TLP.Optimizer()

MOI.set(OPTIMIZER, MOI.Silent(), true)

const CONFIG = MOIT.Config(Float64, atol=1e-6, rtol=1e-6,
    exclude=Any[
        MOI.ConstraintBasisStatus,
        MOI.VariableBasisStatus,
    ]
)

@testset "Direct optimizer" begin

    MOIT.runtests(
        OPTIMIZER, CONFIG,
        exclude=[
            # behaviour to implement: list of model, constraint attributes set
            "test_model_ListOfConstraintAttributesSet",
            "test_model_ModelFilter_AbstractModelAttribute",
            "test_model_ModelFilter_ListOfConstraintIndices",
            "test_model_ModelFilter_ListOfConstraintTypesPresent",
            "test_model_Name",
            "test_objective_set_via_modify",
            # MOI expects to throw when getting duplicate cons / var names
            "test_model_ScalarAffineFunction_ConstraintName",
            "test_model_VariableName",
            "test_model_duplicate_ScalarAffineFunction_ConstraintName",
            "test_model_duplicate_VariableName",
            "test_variable_VariableName",
            # requires get quadratic objective
            "test_objective_get_ObjectiveFunction_ScalarAffineFunction",
            # Tulip not compliant with MOI convention for primal/dual infeasible models
            # See expected behavior at https://jump.dev/MathOptInterface.jl/dev/background/infeasibility_certificates/
            "test_unbounded",
        ]
    )

end

@testset "MOI Bridged" begin
    BRIDGED = MOIB.full_bridge_optimizer(Tulip.Optimizer(), Float64)
    MOI.set(BRIDGED, MOI.Silent(), true)

    MOIT.runtests(
        BRIDGED, CONFIG,
        exclude=[
            # behaviour to implement: list of model, constraint attributes set
            "test_conic_NormInfinityCone_3",
            "test_conic_NormInfinityCone_INFEASIBLE", # should be NO_SOLUTION or INFEASIBLE_POINT
            # ListOfConstraintTypePresent
            "test_conic_NormInfinityCone_VectorAffineFunction",
            "test_conic_NormInfinityCone_VectorOfVariables",
            "test_conic_NormOneCone",
            "test_conic_linear_VectorAffineFunction",
            "test_conic_linear_VectorOfVariables",
            "test_model_delete",
            # List of attributes set
            "test_model_ListOfConstraintAttributesSet",
            "test_model_ModelFilter_AbstractModelAttribute",
            "test_model_ModelFilter_ListOfConstraintIndices",
            "test_model_ModelFilter_ListOfConstraintTypesPresent",
            "test_model_Name",
            "test_objective_set_via_modify",
            # MOI expects to throw when getting duplicate cons / var names
            "test_model_ScalarAffineFunction_ConstraintName",
            "test_model_VariableName",
            "test_model_duplicate_ScalarAffineFunction_ConstraintName",
            "test_model_duplicate_VariableName",
            "test_variable_VariableName",
            # requires get quadratic objective
            "test_objective_get_ObjectiveFunction_ScalarAffineFunction",
            # Tulip not compliant with MOI convention for primal/dual infeasible models
            # See expected behavior at https://jump.dev/MathOptInterface.jl/dev/background/infeasibility_certificates/
            "test_unbounded",
        ],
    )
end

# Run the MOI tests with HSD and MPC algorithms
for ipm in [Tulip.HSD, Tulip.MPC]
    @testset "MOI Linear tests - $ipm" begin
        OPTIMIZER.inner.params.IPM.Factory = Tulip.Factory(ipm)
        MOIT.runtests(OPTIMIZER, CONFIG, include=["linear"])
    end
end

MOIU.@model(ModelData,
        (),
        (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
        (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives),
        (),
        (),
        (MOI.ScalarAffineFunction,),
        (MOI.VectorOfVariables,),
        (MOI.VectorAffineFunction,)
)

@testset "Cached optimizer" begin
    CACHE = MOIU.UniversalFallback(ModelData{Float64}())
    CACHED = MOIU.CachingOptimizer(CACHE, Tulip.Optimizer())
    BRIDGED2 = MOIB.full_bridge_optimizer(CACHED, Float64)
    MOI.set(BRIDGED2, MOI.Silent(), true)

    MOIT.runtests(
        BRIDGED2, CONFIG,
        exclude=[
            # should be NO_SOLUTION or INFEASIBLE_POINT
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
            # Tulip not compliant with MOI convention for primal/dual infeasible models
            # See expected behavior at https://jump.dev/MathOptInterface.jl/dev/background/infeasibility_certificates/
            "test_unbounded",
        ])
end
