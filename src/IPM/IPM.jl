using Printf

"""
    AbstractIPMOptimizer

Abstraction layer for IPM solvers.

An IPM solver implements an interior-point algorithm.
Currently supported:
    * Homogeneous self-dual (HSD)
"""
abstract type AbstractIPMOptimizer{T} end

include("ipmdata.jl")
include("point.jl")
include("residuals.jl")
include("options.jl")


"""
    ipm_optimize!(ipm)

Run the interior-point optimizer of `ipm`.
"""
function ipm_optimize! end

include("HSD/HSD.jl")
include("MPC/MPC.jl")
