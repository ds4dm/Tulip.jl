# Tulip

[![DOI](https://zenodo.org/badge/131298750.svg)](https://zenodo.org/badge/latestdoi/131298750)

 **Documentation** | **Build Status** | **Coverage** |
|:-----------------:|:----------------:|:------------:|
| [![Docs][docs-stable-img]][docs-stable-url] [![Docs-dev][docs-dev-img]][docs-dev-url] | [![Build][build-img]][build-url] [![Build-nightly][build-nightly-img]][build-nightly-url] | [![Codecov][codecov-img]][codecov-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-stable-url]: https://ds4dm.github.io/Tulip.jl/stable
[docs-dev-url]: https://ds4dm.github.io/Tulip.jl/dev/

[build-img]: https://github.com/ds4dm/Tulip.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/ds4dm/Tulip.jl/actions?query=workflow%3ACI
[build-nightly-img]: https://github.com/ds4dm/Tulip.jl/workflows/CI-nightly/badge.svg?branch=master
[build-nightly-url]: https://github.com/ds4dm/Tulip.jl/actions?query=workflow%3A🌙
[codecov-img]: https://codecov.io/github/ds4dm/Tulip.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/ds4dm/Tulip.jl?branch=master


## Overview
Tulip is an open-source interior-point solver for linear optimization, written in pure Julia.
It implements the homogeneous primal-dual interior-point algorithm with multiple centrality corrections, and therefore handles unbounded and infeasible problems.
Tulip’s main feature is that its algorithmic framework is disentangled from linear algebra implementations.
This allows to seamlessly integrate specialized routines for structured problems.

## Installation

Just install like any Julia package

```julia
] add Tulip
```

## Usage

The recommended way of using Tulip is through [JuMP](https://github.com/jump-dev/JuMP.jl) and/or [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) (MOI).

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

The type `Tulip.Optimizer` is parametrized by the model's arithmetic, e.g., `Float64` or `BigFloat`.
This allows to solve problem in higher numerical precision.
See the documentation for more details.

```julia
import MathOptInterface
MOI = MathOptInterface
import Tulip

model = Tulip.Optimizer{Float64}()   # Create a model in Float64 precision
model = Tulip.Optimizer()            # Defaults to the above call
model = Tulip.Optimizer{BigFloat}()  # Create a model in BigFloat precision
```

## Solver parameters

### Setting parameters

When using Tulip through JuMP/MOI, parameters can be set either through MOI's generic `OptimizerAttribute`s, e.g., `MOI.TimeLimitSec` and `MOI.Silent`, or by name.

* Through JuMP
    ```julia
    jump_model = JuMP.Model(Tulip.Optimizer)

    JuMP.set_optimizer_attribute(jump_model, "IPM_IterationsLimit", 200)
    ```

* Through MOI
    ```julia
    moi_model = Tulip.Optimizer{Float64}()

    MOI.set(moi_model, MOI.RawParameter("IPM_IterationsLimit"), 200)
    ```

* Through Tulip's API
    ```julia
    model = Tulip.Model{Float64}()

    Tulip.set_parameter(model, "IPM_IterationsLimit", 200)
    ```

### Parameters description

See the [documentation](https://ds4dm.github.io/Tulip.jl/stable/reference/options/).

## Citing `Tulip.jl`

If you use Tulip in your work, we kindly ask that you cite the following reference.
The PDF is freely available [here](https://www.gerad.ca/fr/papers/G-2019-36/view), and serves as a user manual for advanced users.

```
@TechReport{Tulip.jl,
    title = {{Tulip}.jl: an open-source interior-point linear optimization
    solver with abstract linear algebra},
    url = {https://www.gerad.ca/fr/papers/G-2019-36},
    Journal = {Les Cahiers du Gerad},
    Author = {Anjos, Miguel F. and Lodi, Andrea and Tanneau, Mathieu},
    year = {2019}
}
```
