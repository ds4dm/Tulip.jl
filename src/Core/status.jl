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
    # OK statuses
    Trm_Optimal,
    Trm_PrimalInfeasible,
    Trm_DualInfeasible,
    Trm_PrimalDualInfeasible,
    # Limits
    Trm_IterationLimit,
    Trm_TimeLimit,
    # Errors
    Trm_MemoryLimit,
    Trm_NumericalProblem
)

"""
    SolutionStatus

Solution Status code
- `Sln_Unknown`: Unknown status
- `Sln_Optimal`: The current solution is optimal.
- `Sln_FeasiblePoint`: The current solution is feasible.
- `Sln_InfeasiblePoint`: The current solution is not feasible.
- `Sln_InfeasibilityCertificate`: The current solution is a certificate of 
    infeasibility. The primal solution is a certificate of dual infeasibility, 
    while the dual solution is a certificate of primal infeasibility.
"""
@enum(SolutionStatus,
    Sln_Unknown,
    Sln_Optimal,
    Sln_FeasiblePoint,
    Sln_InfeasiblePoint,
    Sln_InfeasibilityCertificate
)