const _KRYLOV_SPD = Union{
    Krylov.CgSolver,
    Krylov.CrSolver,
}

const _KRYLOV_SID = Union{
    Krylov.MinresSolver,
    Krylov.MinresQlpSolver,
    Krylov.SymmlqSolver
}

const _KRYLOV_SQD = Union{
    Krylov.TricgSolver,
    Krylov.TrimrSolver,
}

const _KRYLOV_LN = Union{
    Krylov.LnlqSolver,
    Krylov.CraigSolver,
    Krylov.CraigmrSolver,
}

const _KRYLOV_LS = Union{
    Krylov.LslqSolver,
    Krylov.LsqrSolver,
    Krylov.LsmrSolver,
}

# Helper functions
for (KS, fun) in [
    (Krylov.CgSolver,Krylov.cg!)
    (Krylov.CrSolver,Krylov.cr!)
    (Krylov.MinresSolver,Krylov.minres!)
    (Krylov.MinresQlpSolver,Krylov.minres_qlp!)
    (Krylov.SymmlqSolver,Krylov.symmlq!)
    (Krylov.TricgSolver,Krylov.tricg!)
    (Krylov.TrimrSolver,Krylov.trimr!)
    (Krylov.LnlqSolver,Krylov.lnlq!)
    (Krylov.CraigSolver,Krylov.craig!)
    (Krylov.CraigmrSolver,Krylov.craigmr!)
    (Krylov.LslqSolver,Krylov.lslq!)
    (Krylov.LsqrSolver,Krylov.lsqr!)
    (Krylov.LsmrSolver,Krylov.lsmr!)
]
    @eval begin
        @inline _krylov!(solver::$KS, args...; kwargs...) = $(fun)(solver, args...; kwargs...)
    end
end
