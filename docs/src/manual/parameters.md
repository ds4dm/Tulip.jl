# Parameters

## Accessing and modifying parameters

### Directly

Parameters are stored in the `env` attribute of a `Model` object.
To access/modify a parameter, simply access the corresponding field as follows
```julia
model = Tulip.Model{Float64}()

model.env.barrier_tol_pfeas = 1e-10  # primal feasibility tolerance
model.env.barrier_tol_dfeas = 1e-10  # dual feasibility tolerance
```

See [List of parameters](@ref) below for a more exhaustive list of available parameters.

### Through JuMP / MathOptInterface

Some parameters can be accessed/modified through Tulip's `MathoptInterface` API.
Currently, only `MOI.TimeLimitSec`, `MOI.NumberOfThreads`, `MOI.Silent` are supported.

Managing other parameters through `MOI.RawParameter` and `JuMP.set_optimizer_attribute` is not supported at this time.
Therefore, one needs to access the `Tulip.Model` object and set parameters directly as follows
```julia
moi_model = Tulip.Optimizer{Float64}()
tlp_model = moi_model.inner

# Set some parameters
tlp_model.env.barrier_tol_pfeas = 1e-10  # primal feasibility tolerance
tlp_model.env.barrier_tol_dfeas = 1e-10  # dual feasibility tolerance
```

## List of parameters

### Tolerances

Numerical tolerances for the interior-point algorithm.

| Parameter | Description | Default |
|:----------|:------------|--------:|
| `barrier_tol_pfeas` | primal feasibility tolerance | ``\sqrt{\epsilon}``
| `barrier_tol_dfeas` | dual feasibility tolerance | ``\sqrt{\epsilon}``
| `barrier_tol_conv` | optimality tolerance | ``\sqrt{\epsilon}``
| `barrier_tol_infeas` | infeasibility tolerance | ``\sqrt{\epsilon}``

where ``\epsilon`` is the numerical precision.

, which is given by `eps(Tv)`.
For `Tv=Float64`, this yields ``\epsilon \simeq 10^{-16}`` and thus tolerances of about ``10^{-8}``.


### Stopping criterion

| Parameter | Description | Default |
|:----------|:------------|--------:|
| `barrier_iter_max` | maximum number of barrier iterations | `100` |
| `time_limit` | time limit, in seconds | `Inf` |


### Others

| Parameter | Description | Default |
|:----------|:------------|--------:|
| `verbose` | verbosity flag | `0` |
| `threads` | maximum number of threads | `1` |