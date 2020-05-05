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
    
    # Algorithmic parameters
    BarrierAlgorithm::Int = 1  # TODO: docs
    BarrierCorrectionLimit::Int = 5  # Maximum number of centrality corrections
    BarrierStepDampFactor::Tv = Tv(9_995 // 10_000)  # Damp step size by this much
    BarrierGammaMin::Tv = Tv(1 // 10)
    BarrierCentralityOutlierThreshold::Tv = Tv(1 // 10)  # Relative threshold for centrality outliers

    # Linear algebra
    MatrixOptions::TLA.MatrixOptions = TLA.MatrixOptions(SparseMatrixCSC)
    KKTOptions::KKT.SolverOptions = KKT.default_options(Tv)
    BarrierPRegMin::Tv = sqrt(eps(Tv))  # primal
    BarrierDRegMin::Tv = sqrt(eps(Tv))  # dual
    # TODO: iterative refinement

    # I/O
    OutputLevel::Int = 0
    LogLevel::Int = 0  # TODO: have different LogLevel for, e.g., presolve, optimize!, etc...

    # Others
    Threads::Int = 1
    Presolve::Int = 1
end