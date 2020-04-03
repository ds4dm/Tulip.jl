```@meta
CurrentModule = Tulip
```

# Parameters

Parameters can be queried and set through the [`get_parameter`](@ref) and [`set_parameter`](@ref) functions.


In all that follows, ``\epsilon`` refers to the numerical precision, which is given by `eps(Tv)` where `Tv` is the arithmetic of the current model.
For instance, in double-precision arithmetic, i.e., `Tv=Float64`, we have ``\epsilon \simeq 10^{-16}``.

## Tolerances

Numerical tolerances for the interior-point algorithm.

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `BarrierTolerancePFeas` | Primal feasibility tolerance | ``\sqrt{\epsilon}``
| `BarrierToleranceDFeas` | Dual feasibility tolerance | ``\sqrt{\epsilon}``
| `BarrierToleranceRGap`  | Relative optimality gap tolerance | ``\sqrt{\epsilon}``
| `BarrierToleranceIFeas` | Infeasibility tolerance | ``\sqrt{\epsilon}``

## Algorithmic parameters

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `BarrierAlgorithm` | Interior-point algorithm | `1` |
| `BarrierCorrectionLimit` | Maximum number of additional centrality corrections | `5` |
| `BarrierStepDampFactor` | Step | `0.9995` |
| `BarrierGammaMin` | Minimum value of ``\gamma`` for computing correctors | `0.1`
| `BarrierCentralityOutlierThreshold` | Relative threshold for computing extra centrality corrections | `0.1`
| `BarrierPRegMin` | Minimum value of primal regularization | ``\sqrt{\epsilon}`` |
| `BarrierDregMin` | Minimum value of dual regularization | ``\sqrt{\epsilon}``

## Stopping criterion

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `BarrierIterationsLimit` | Maximum number of barrier iterations | `100` |
| `TimeLimit` | Time limit, in seconds | `Inf` |

## Linear solvers

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `LinearSolverBackend` | Backend for solving linear systems | [`Tulip.TLPLinearAlgebra.DefaultBackend`](@ref) |
| `LinearSolverSystem` | Linear system solved | [`Tulip.TLPLinearAlgebra.DefaultSystem`](@ref) |

## Others

| Parameter | Description | Default |
|:----------|:------------|:--------|
| `OutputLevel` | Controls the solver's output | `0` |
| `Threads` | Maximum number of threads | `1` |