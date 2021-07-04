const KrylovSPD = Union{
    Krylov.CgSolver,
    Krylov.CrSolver,
    Krylov.MinresSolver,
    Krylov.MinresQlpSolver,
    Krylov.SymmlqSolver
}

const KrylovSID = Union{
    Krylov.MinresSolver,
    Krylov.MinresQlpSolver,
    Krylov.SymmlqSolver
}

# Helper functions
# TODO: use macros to generate those automatically
@inline _krylov!(ksolver::CgSolver, A, b; kwargs...) = Krylov.cg!(ksolver, A, b; kwargs...)
@inline _krylov!(ksolver::CrSolver, A, b; kwargs...) = Krylov.cr!(ksolver, A, b; kwargs...)
@inline _krylov!(ksolver::MinresSolver, A, b; kwargs...) = Krylov.minres!(ksolver, A, b; kwargs...)
@inline _krylov!(ksolver::MinresQlpSolver, A, b; kwargs...) = Krylov.minres_qlp!(ksolver, A, b; kwargs...)
@inline _krylov!(ksolver::SymmlqSolver, A, b; kwargs...) = Krylov.symmlq!(ksolver, A, b; kwargs...)
