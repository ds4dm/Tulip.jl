```@meta
CurrentModule = Tulip.TLPLinearAlgebra
```

### AbstractLinearSolver

This is the base type from which all implementations should derive.

```@docs
AbstractLinearSolver
```

Custom linear solvers should inherit from the `AbstractLinearSolver` class,
and extend the following two functions:

```@docs
update_linear_solver!
```

```@docs
solve_augmented_system!
```

### DenseLinearSolver

```@docs
DenseLinearSolver
```

```@docs
update_linear_solver!(::DenseLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
update_linear_solver!(::DenseLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:BlasReal}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::DenseLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::DenseLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:BlasReal}
```

### SparseIndefLinearSolver

```@docs
SparseIndefLinearSolver
```

```@docs
update_linear_solver!(::SparseIndefLinearSolver,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::SparseIndefLinearSolver, ::Vector{Float64}, ::Vector{Float64})
```



### SparsePosDefLinearSolver

```@docs
SparsePosDefLinearSolver
```

```@docs
update_linear_solver!(::SparsePosDefLinearSolver,::AbstractVector{Float64},::AbstractVector{Float64},::AbstractVector{Float64})
```

```@docs
solve_augmented_system!(::Vector{Float64},::Vector{Float64},::SparsePosDefLinearSolver, ::Vector{Float64}, ::Vector{Float64})
```

### LDLFLinearSolver

```@docs
LDLFLinearSolver
```

```@docs
update_linear_solver!(::LDLFLinearSolver{Tv},::AbstractVector{Tv},::AbstractVector{Tv},::AbstractVector{Tv}) where{Tv<:Real}
```

```@docs
solve_augmented_system!(::Vector{Tv},::Vector{Tv},::LDLFLinearSolver{Tv}, ::Vector{Tv}, ::Vector{Tv}) where{Tv<:Real}
```