mutable struct Solution{T}
    m::Int
    n::Int

    primal_status::SolutionStatus
    dual_status::SolutionStatus
    is_primal_ray::Bool
    is_dual_ray::Bool

    z_primal::T
    z_dual::T

    x::Vector{T}
    Ax::Vector{T}
    y_lower::Vector{T}
    y_upper::Vector{T}
    s_lower::Vector{T}
    s_upper::Vector{T}

    Solution{T}(m, n) where{T} = new{T}(
        m, n, Sln_Unknown, Sln_Unknown, false, false,
        zero(T), zero(T),
        zeros(T, n), zeros(T, m),
        zeros(T, m), zeros(T, m),
        zeros(T, n), zeros(T, n)
    )
end

import Base.resize!


function Base.resize!(sol::Solution, m::Int, n::Int)
    m >= 0 || throw(ArgumentError("m must be >= 0"))
    n >= 0 || throw(ArgumentError("n must be >= 0"))
    
    sol.m = m
    sol.n = n

    resize!(sol.x, n)
    resize!(sol.Ax, m)

    resize!(sol.y_lower, m)
    resize!(sol.y_upper, m)
    resize!(sol.s_lower, n)
    resize!(sol.s_upper, n)

    return sol
end