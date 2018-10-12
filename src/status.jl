"""
    TerminationStatus

- `Success`: No error occured
- `PrimalInfeasibleNoResult`: Problem is proved to be primal infeasible,
    but no result (e.g. certificate of infeasibility) is available.
- `DualInfeasibleNoResult`: Problem is proved to be primal infeasible,
but no result (e.g. certificate of infeasibility) is available.
- `IterationLimit`: Maximum number of iterations reached.
- `TimeLimit`: Time limit reached.
- `MemoryLimit`: Memory limit reached.
- `NumericalProblem`: Numerical problem encountered, e.g. failure of the
    Cholesky decomposition.
"""
@enum(TerminationStatus,
    Trm_Unknown,
    Trm_Success,
    Trm_PrimalInfeasibleNoResult,
    Trm_DualInfeasibleNoResult,
    Trm_IterationLimit,
    Trm_TimeLimit,
    Trm_MemoryLimit,
    Trm_NumericalProblem
)

"""
    SolutionStatus

Solution Status code
- `Unknown`: Unknown status
- `Optimal`: A proven optimal solution is available.
- `PrimalFeasible`: A feasible primal solution is available.
- `DualFeasible`: A feasible dual solution is available.
- `PrimalDualFeasible`: A feasible primal-dual solution pair is available.
- `PrimalInfeasible`: Problem was proved to be primal infeasible, and an infeasibility
    certificate is available.
- `DualInfeasible`: Problem was proved to be dual infeasible, and an infeasibility
    certificate is available.
- `PrimalDualInfeasible`: Both primal and dual are infeasible, and both infeasiblity
    certificate are available.
"""
@enum(SolutionStatus,
    Sln_Unknown,
    Sln_Optimal,
    Sln_PrimalFeasible,
    Sln_DualFeasible,
    Sln_PrimalDualFeasible,
    Sln_PrimalInfeasible,
    Sln_DualInfeasible,
    Sln_PrimalDualInfeasible
)

"""
    ModelStatus

Model status codes.
- `Mdl_Unknown`: Status of the model is unknown.
- `Mdl_StandardForm`: Model is in standard form.
"""
@enum(ModelStatus,
    Mdl_Unknown,
    Mdl_StandardForm
)