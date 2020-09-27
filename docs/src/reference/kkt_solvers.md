```@meta
CurrentModule = Tulip.KKT
```

## Overview

### AbstractKKTSolver

This is the base type from which all implementations should derive.

```@docs
AbstractKKTSolver
```

Custom linear solvers should inherit from the `AbstractKKTSolver` class,
and extend the following two functions:

```@docs
update!
```

```@docs
solve!
```

### Choosing between linear solvers

```@docs
SolverOptions
```

```@meta
CurrentModule = Tulip
```

```@docs
TLA.MatrixOptions
```

```@meta
CurrentModule = Tulip.KKT
```

## Dense/LAPACK

```@docs
DenseSPD
```

## CHOLMOD

```@docs
CholmodSolver
```

```@docs
CholmodSQD
```

```@docs
CholmodSPD
```

## LDLFactorizations

```@docs
LDLFactSQD
```

## Krylov

!!! warning
    Iterative methods are still an experimental feature.
    Some numerical and performance issues should be expected.


```@docs
KrylovSPD
```

```@docs
KrylovSID
```

```@docs
KrylovSQD
```