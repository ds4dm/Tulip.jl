# Tulip

[![DOI](https://zenodo.org/badge/131298750.svg)](https://zenodo.org/badge/latestdoi/131298750)
[![](https://github.com/ds4dm/Tulip.jl/workflows/CI/badge.svg?branch=master)](https://github.com/ds4dm/Tulip.jl/actions?query=workflow%3ACI)
[![](https://codecov.io/github/ds4dm/Tulip.jl/coverage.svg?branch=master)](https://codecov.io/github/ds4dm/Tulip.jl?branch=master)

[Tulip](https://github.com/ds4dm/Tulip.jl) is an open-source interior-point solver for linear optimization, written in pure Julia.
It implements the homogeneous primal-dual interior-point algorithm with multiple centrality corrections, and therefore handles unbounded and infeasible problems.
Tulip’s main feature is that its algorithmic framework is disentangled from linear algebra implementations.
This allows to seamlessly integrate specialized routines for structured problems.

## License

Tulip is licensed under the [MPL 2.0 license](https://github.com/ds4dm/Tulip.jl/blob/master/LICENSE.md).

## Installation

Install Tulip using the Julia package manager:

```julia
import Pkg
Pkg.add("Tulip")
```

## Usage

The recommended way of using Tulip is through [JuMP](https://github.com/jump-dev/JuMP.jl) or [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) (MOI).

The low-level interface is still under development and is likely change in the future.
The MOI interface is more stable.

### Using with JuMP

Tulip follows the syntax convention `PackageName.Optimizer`:

```julia
using JuMP
import Tulip
model = Model(Tulip.Optimizer)
```

Linear objectives, linear constraints and lower/upper bounds on variables are supported.

### Using with MOI

The type `Tulip.Optimizer` is parametrized by the model's arithmetic, for example, `Float64` or `BigFloat`.
This allows to solve problem in higher numerical precision.
See the documentation for more details.

```julia
import MathOptInterface as MOI
import Tulip
model = Tulip.Optimizer{Float64}()   # Create a model in Float64 precision
model = Tulip.Optimizer()            # Defaults to the above call
model = Tulip.Optimizer{BigFloat}()  # Create a model in BigFloat precision
```

## Solver parameters

See the [documentation](https://ds4dm.github.io/Tulip.jl/stable/reference/options/) for a full list of parameters.

To set parameters in JuMP, use:
```julia
using JuMP, Tulip
model = Model(Tulip.Optimizer)
set_attribute(model, "IPM_IterationsLimit", 200)
```

To set parameters in MathOptInterface, use:
```julia
using Tulip
import MathOptInterface as MOI
model = Tulip.Optimizer{Float64}()
MOI.set(model, MOI.RawOptimizerAttribute("IPM_IterationsLimit"), 200)
```

To set parameters in the Tulip API, use:
```julia
using Tulip
model = Tulip.Model{Float64}()
Tulip.set_parameter(model, "IPM_IterationsLimit", 200)
```

## Command-line executable

See [app building instructions](https://github.com/ds4dm/Tulip.jl/blob/master/app/README.md).

## Citing `Tulip.jl`

If you use Tulip in your work, we kindly ask that you cite the following [reference](https://doi.org/10.1007/s12532-020-00200-8) (preprint available [here](https://arxiv.org/abs/2006.08814)).

```
@Article{Tulip.jl,
  author   = {Tanneau, Mathieu and Anjos, Miguel F. and Lodi, Andrea},
  journal  = {Mathematical Programming Computation},
  title    = {Design and implementation of a modular interior-point solver for linear optimization},
  year     = {2021},
  issn     = {1867-2957},
  month    = feb,
  doi      = {10.1007/s12532-020-00200-8},
  language = {en},
  url      = {https://doi.org/10.1007/s12532-020-00200-8},
  urldate  = {2021-03-07},
}
```
