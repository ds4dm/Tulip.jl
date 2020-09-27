# TODO: enum for some parameters?
# Or something easier to expand?
"""
    Parameters{T}

"""
Base.@kwdef mutable struct Parameters{T}
    
    # User limits
    BarrierIterationsLimit::Int = 100
    TimeLimit::Float64 = Inf

    # Numerical tolerances
    BarrierTolerancePFeas::T = sqrt(eps(T))  # primal feasibility
    BarrierToleranceDFeas::T = sqrt(eps(T))  # dual feasibility
    BarrierToleranceRGap::T = sqrt(eps(T))  # optimality
    BarrierToleranceIFeas::T = sqrt(eps(T))  # infeasibility
    
    # Algorithmic parameters
    BarrierAlgorithm::Int = 1  # TODO: docs
    BarrierCorrectionLimit::Int = 5  # Maximum number of centrality corrections
    BarrierStepDampFactor::T = T(9_995 // 10_000)  # Damp step size by this much
    BarrierGammaMin::T = T(1 // 10)
    BarrierCentralityOutlierThreshold::T = T(1 // 10)  # Relative threshold for centrality outliers

    # Linear algebra
    MatrixOptions::TLA.MatrixOptions = TLA.MatrixOptions(SparseMatrixCSC)
    KKTOptions::KKT.SolverOptions = KKT.default_options(T)
    BarrierPRegMin::T = sqrt(eps(T))  # primal
    BarrierDRegMin::T = sqrt(eps(T))  # dual
    # TODO: iterative refinement

    # I/O
    OutputLevel::Int = 0
    LogLevel::Int = 0  # TODO: have different LogLevel for, e.g., presolve, optimize!, etc...

    # Others
    Threads::Int = 1
    Presolve::Int = 1
end