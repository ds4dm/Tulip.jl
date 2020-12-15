Base.@kwdef mutable struct IPMOptions{T}

    OutputLevel::Int = 0

    # User limits
    IterationsLimit::Int = 100
    TimeLimit::Float64 = Inf

    # Numerical tolerances
    TolerancePFeas::T = sqrt(eps(T))  # primal feasibility
    ToleranceDFeas::T = sqrt(eps(T))  # dual feasibility
    ToleranceRGap::T = sqrt(eps(T))  # optimality
    ToleranceIFeas::T = sqrt(eps(T))  # infeasibility

    # Algorithmic parameters
    CorrectionLimit::Int = 3  # Maximum number of centrality corrections
    StepDampFactor::T = T(9_995 // 10_000)  # Damp step size by this much
    GammaMin::T = T(1 // 10)
    CentralityOutlierThreshold::T = T(1 // 10)  # Relative threshold for centrality outliers

    PRegMin::T = sqrt(eps(T))  # primal
    DRegMin::T = sqrt(eps(T))  # dual

    Factory::Factory{<:AbstractIPMOptimizer} = Factory(HSD)
end
