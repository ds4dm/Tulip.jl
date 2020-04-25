# TODO: enum for some parameters?
# Or something easier to expand?
"""
    Parameters{Tv}

"""
Base.@kwdef mutable struct Parameters{Tv}
    
    # User limits
    BarrierIterationsLimit::Int = 100
    TimeLimit::Float64 = Inf

    # Numerical tolerances
    BarrierTolerancePFeas::Tv = sqrt(eps(Tv))  # primal feasibility
    BarrierToleranceDFeas::Tv = sqrt(eps(Tv))  # dual feasibility
    BarrierToleranceRGap::Tv = sqrt(eps(Tv))  # optimality
    BarrierToleranceIFeas::Tv = sqrt(eps(Tv))  # infeasibility

    # Regularizations
    BarrierPRegMin::Tv = sqrt(eps(Tv))  # primal
    BarrierDregMin::Tv = sqrt(eps(Tv))  # dual
    
    # Algorithmic parameters
    BarrierAlgorithm::Int = 1  # TODO: docs
    BarrierCorrectionLimit::Int = 5  # Maximum number of centrality corrections
    BarrierStepDampFactor::Tv = Tv(9_995 // 10_000)  # Damp step size by this much
    BarrierGammaMin::Tv = Tv(1 // 10)
    BarrierCentralityOutlierThreshold::Tv = Tv(1 // 10)  # Relative threshold for centrality outliers

    # Linear algebra
    MatrixType::Type{<:AbstractMatrix} = SparseMatrixCSC  # TODO: make this more flexible
    #=
        TODO: make it possible to use custom linear solver,
            without having to change the source code.
        KKTBackend:
        * -1: Default
        *  0: LAPACK (dense)
        *  1: Julia's generic dense Cholesky
        *  2: CHOLMOD (sparse)
        *  3: LDLFact (sparse)
    
        KKTSystem:
        * -1: Default
        *  0: Augmented system
        *  1: Normal equations
    =#
    KKTBackend::Int = -1
    KKTSystem::Int = -1
    LinearSolverBackend::TLA.LSBackend = TLA.DefaultBackend()
    LinearSolverSystem::TLA.LinearSystem = TLA.DefaultSystem()
    # TODO: put regularizations here (?)
    # TODO: iterative refinement

    # I/O
    OutputLevel::Int = 0
    LogLevel::Int = 0  # TODO: have different LogLevel for, e.g., presolve, optimize!, etc...

    # Others
    Threads::Int = 1
    Presolve::Int = 1
end