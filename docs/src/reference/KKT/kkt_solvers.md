```@meta
CurrentModule = Tulip.KKT
```
# KKT solvers

## Overview

Different KKT solvers can be selected via a combination of _backend_ and _KKT system_.

Here is a short overview of which backend supports which systems & arithmetic

| Backend                       | Systems                    | Arithmetic   | Method         |
|:------------------------------|:---------------------------|:-------------|:---------------|
| [`TlpCholmod.Backend`](@ref)  | [`K2`](@ref), [`K1`](@ref) | `Float64`    | LDLᵀ, Cholesky
| [`TlpDense.Backend`](@ref)    | [`K1`](@ref)               | `Any`        | Cholesky
| [`TlpLDLFact.Backend`](@ref)  | [`K2`](@ref)               | `Any`        | LDLᵀ
| [`TlpKrylov.Backend`](@ref)   | [`K2`](@ref), [`K1`](@ref) | `Any`        | Krylov

### AbstractKKTSolver

This is the base type from which all implementations should derive.

```@docs
AbstractKKTSolver
```

Custom linear solvers should (preferably) inherit from the `AbstractKKTSolver` class,
and extend the following functions:

```@docs
setup
```

```@docs
update!
```

```@docs
solve!
```

### Choosing between linear solvers

```@docs
KKTOptions
```

## `TlpCholmod`

```@docs
TlpCholmod.Backend
```

```@docs
TlpCholmod.CholmodSolver
```

## `TlpDense`

```@docs
TlpDense.Backend
```

```@docs
TlpDense.DenseSolver
```

## `TlpLDLFact`

```@docs
TlpLDLFactorizations.Backend
```

```@docs
TlpLDLFactorizations.LDLFactSolver
```

## `TlpKrylov`

!!! warning
    Iterative methods are still an experimental feature.
    Some numerical and performance issues should be expected.


```@docs
TlpKrylov.Backend
```

```@docs
TlpKrylov.AbstractKrylovSolver
```

```@docs
TlpKrylov.SPDSolver
```

```@docs
TlpKrylov.SIDSolver
```

```@docs
TlpKrylov.SQDSolver
```