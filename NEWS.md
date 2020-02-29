# Tulip release notes

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