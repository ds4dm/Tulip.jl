Base.@kwdef mutable struct IPMOptions{T}

    OutputLevel::Int = 0

    # User limits
    BarrierIterationsLimit::Int = 100
    TimeLimit::Float64 = Inf

    # Numerical tolerances
    BarrierTolerancePFeas::T = sqrt(eps(T))  # primal feasibility
    BarrierToleranceDFeas::T = sqrt(eps(T))  # dual feasibility
    BarrierToleranceRGap::T = sqrt(eps(T))  # optimality
    BarrierToleranceIFeas::T = sqrt(eps(T))  # infeasibility

    # Algorithmic parameters
    BarrierCorrectionLimit::Int = 5  # Maximum number of centrality corrections
    BarrierStepDampFactor::T = T(9_995 // 10_000)  # Damp step size by this much
    BarrierGammaMin::T = T(1 // 10)
    BarrierCentralityOutlierThreshold::T = T(1 // 10)  # Relative threshold for centrality outliers

    BarrierPRegMin::T = sqrt(eps(T))  # primal
    BarrierDRegMin::T = sqrt(eps(T))  # dual

    Factory::Factory{<:AbstractIPMOptimizer} = Factory(HSD)
end