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
update!
```

```@docs
solve!
```

### KKTSolver_Dense

```@docs
KKTSolver_Dense
```

```@docs
update!(::KKTSolver_Dense{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
update!(::KKTSolver_Dense{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:BlasReal}
```

```@docs
solve!(::Vector{Tv},::Vector{Tv},::KKTSolver_Dense{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
solve!(::Vector{Tv},::Vector{Tv},::KKTSolver_Dense{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:BlasReal}
```

### KKTSolver_CholmodQD

```@docs
KKTSolver_CholmodQD
```

```@docs
update!(::KKTSolver_CholmodQD,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve!(::Vector{Float64},::Vector{Float64},::KKTSolver_CholmodQD, ::Vector{Float64}, ::Vector{Float64})
```



### KKTSolver_CholmodPD

```@docs
KKTSolver_CholmodPD
```

```@docs
update!(::KKTSolver_CholmodPD,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve!(::Vector{Float64},::Vector{Float64},::KKTSolver_CholmodPD, ::Vector{Float64}, ::Vector{Float64})
```

### KKTSolver_LDLFact

```@docs
KKTSolver_LDLFact
```

```@docs
update!(::KKTSolver_LDLFact{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
```

```@docs
solve!(::Vector{Tv},::Vector{Tv},::KKTSolver_LDLFact{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
```