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
    Success,
    PrimalInfeasibleNoResult,
    DualInfeasibleNoResult,
    IterationLimit,
    TimeLimit,
    MemoryLimit,
    NumericalProblem
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
    Unknown,
    Optimal,
    PrimalFeasible,
    DualFeasible,
    PrimalDualFeasible,
    PrimalInfeasible,
    DualInfeasible,
    PrimalDualInfeasible
)

