```@meta
CurrentModule = Tulip.KKT
```

### AbstractKKTSolver

This is the base type from which all implementations should derive.

```@docs
AbstractKKTSolver
```

Custom linear solvers should inherit from the `AbstractKKTSolver` class,
and extend the following two functions:

```@docs
update_linear_solver!
```

```@docs
solve_augmented_system!
```

### KKTSolver_Dense

```@docs
KKTSolver_Dense
```

```@docs
update_linear_solver!(::KKTSolver_Dense{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
update_linear_solver!(::KKTSolver_Dense{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:BlasReal}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::KKTSolver_Dense{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::KKTSolver_Dense{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:BlasReal}
```

### KKTSolver_CholmodQD

```@docs
KKTSolver_CholmodQD
```

```@docs
update_linear_solver!(::KKTSolver_CholmodQD,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::KKTSolver_CholmodQD, ::Vector{Float64}, ::Vector{Float64})
```



### KKTSolver_CholmodPD

```@docs
KKTSolver_CholmodPD
```

```@docs
update_linear_solver!(::KKTSolver_CholmodPD,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::KKTSolver_CholmodPD, ::Vector{Float64}, ::Vector{Float64})
```

### KKTSolver_LDLFact

```@docs
KKTSolver_LDLFact
```

```@docs
update_linear_solver!(::KKTSolver_LDLFact{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::KKTSolver_LDLFact{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
```