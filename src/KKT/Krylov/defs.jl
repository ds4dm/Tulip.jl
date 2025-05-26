if isdefined(Krylov, :CgWorkspace)
    const _Krylov = Krylov
else
    @eval module _Krylov
    using Krylov
    const MinresWorkspace = Krylov.MinresSolver
    const CgWorkspace = Krylov.CgSolver
    const CrWorkspace = Krylov.CrSolver
    const SymmlqWorkspace = Krylov.SymmlqSolver
    const MinresQlpWorkspace = Krylov.MinresQlpSolver
    const TricgWorkspace = Krylov.TricgSolver
    const TrimrWorkspace = Krylov.TrimrSolver
    const LslqWorkspace = Krylov.LslqSolver
    const LsqrWorkspace = Krylov.LsqrSolver
    const LsmrWorkspace = Krylov.LsmrSolver
    const LnlqWorkspace = Krylov.LnlqSolver
    const CraigWorkspace = Krylov.CraigSolver
    const CraigmrWorkspace = Krylov.CraigmrSolver
    end
end

const _KRYLOV_SPD = Union{
    _Krylov.CgWorkspace,
    _Krylov.CrWorkspace,
}

const _KRYLOV_SID = Union{
    _Krylov.MinresWorkspace,
    _Krylov.MinresQlpWorkspace,
    _Krylov.SymmlqWorkspace
}

const _KRYLOV_SQD = Union{
    _Krylov.TricgWorkspace,
    _Krylov.TrimrWorkspace,
}

const _KRYLOV_LN = Union{
    _Krylov.LnlqWorkspace,
    _Krylov.CraigWorkspace,
    _Krylov.CraigmrWorkspace,
}

const _KRYLOV_LS = Union{
    _Krylov.LslqWorkspace,
    _Krylov.LsqrWorkspace,
    _Krylov.LsmrWorkspace,
}

# Helper functions
for (KS, fun) in [
    (_Krylov.CgWorkspace,Krylov.cg!)
    (_Krylov.CrWorkspace,Krylov.cr!)
    (_Krylov.MinresWorkspace,Krylov.minres!)
    (_Krylov.MinresQlpWorkspace,Krylov.minres_qlp!)
    (_Krylov.SymmlqWorkspace,Krylov.symmlq!)
    (_Krylov.TricgWorkspace,Krylov.tricg!)
    (_Krylov.TrimrWorkspace,Krylov.trimr!)
    (_Krylov.LnlqWorkspace,Krylov.lnlq!)
    (_Krylov.CraigWorkspace,Krylov.craig!)
    (_Krylov.CraigmrWorkspace,Krylov.craigmr!)
    (_Krylov.LslqWorkspace,Krylov.lslq!)
    (_Krylov.LsqrWorkspace,Krylov.lsqr!)
    (_Krylov.LsmrWorkspace,Krylov.lsmr!)
    ]
    @eval begin
        @inline _krylov!(solver::$KS, args...; kwargs...) = $(fun)(solver, args...; kwargs...)
    end
end
