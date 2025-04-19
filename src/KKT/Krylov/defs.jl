const _KRYLOV_SPD = Union{
    Krylov.CgWorkspace,
    Krylov.CrWorkspace,
    Krylov.CarWorkspace,
}

const _KRYLOV_SID = Union{
    Krylov.MinresWorkspace,
    Krylov.MinaresWorkspace,
    Krylov.MinresQlpWorkspace,
    Krylov.SymmlqWorkspace
}

const _KRYLOV_SQD = Union{
    Krylov.TricgWorkspace,
    Krylov.TrimrWorkspace,
}

const _KRYLOV_LN = Union{
    Krylov.LnlqWorkspace,
    Krylov.CraigWorkspace,
    Krylov.CraigmrWorkspace,
}

const _KRYLOV_LS = Union{
    Krylov.LslqWorkspace,
    Krylov.LsqrWorkspace,
    Krylov.LsmrWorkspace,
}
