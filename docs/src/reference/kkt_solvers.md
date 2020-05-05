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
Dense_SymPosDef
```

## CHOLMOD

```@docs
CholmodSolver
```

```@docs
Cholmod_SymQuasDef
```

```@docs
Cholmod_SymPosDef
```

## LDLFactorizations

```@docs
LDLFact_SymQuasDef
```