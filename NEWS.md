# Tulip release notes

# v0.7.1 (January 8, 2021)

## Bug fixes
* Fix a bug in `add_variable!` and the handling of zero coefficients(#77, #79)
* Fix some type instabilities (#71, #76)
* Fix docs URL in README (#70)

## New features
* Tulip's version is exposed via `Tulip.version()` (#81)
* Multiple centrality corrections in MPC algorithm (#75)
* Support reading `.gz` and `.bz2` files (#72)

## Others
* Move CI to GitHub actions (#73)
* Use new convention for test-specific dependencies (#80)
* Bump deps (#71,#74)

# v0.5.1 (August 15, 2020)
* Fix URLs following migration to jump.dev (#51)
* Fix bug in constraint/variable deletion (#52, #53)

# v0.5.0
* Support the MOI attribute `Name` (#47)
* Simplify the user interface for choosing between different linear solvers (#48)

# v0.4.0 (April 25, 2020)
* Re-write data structure and interface (#44)
* Add presolve module (#45)
* Move `UnitBlockAngular` code to a separate package (#46)

# v0.3.0 (February 29, 2020)

* Improved documentation for parameter management (#43)
* More flexible management of linear solvers (#37, #40)
    * Introduce two new parameters to choose linear solver
        * `ls_backend` specifies which backend is used to solve linear systems
        * `ls_system` specifies which linear system (augmented system or normal equations) is solved
    * Generic tests for custom linear solvers
    * Performance improvements when using CHOLMOD
* Improved MPS reader (#36)

# v0.2.0 (December 26, 2019)

* Switch type unions to constant in the MOI wrapper (#32)
* Free MPS format reader (#33).
Some `.mps` files in fixed MPS format may no longer be readable, due to, e.g., 
    * spaces in constraint/variable names
    * empty name field in RHS section
* Re-write of the linear algebra layer (#34)
* Integration of [LDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LDLFactorizations.jl) (#35) for solving problems in arbitrary precision.

# v0.1.1 (October 25, 2019)

* Support for MOI v0.9.5 (#31)
* Catch exceptions during solving (#29)
    * Numerical trouble during factorization
    * Iteration and memory limits
    * User interruptions
* Create a `CITATION.bib` file (#26)
* Describe problem formulations in the docs (#26)