mutable struct Solution{Tv}
    m::Int
    n::Int

    primal_status::SolutionStatus
    dual_status::SolutionStatus
    is_primal_ray::Bool
    is_dual_ray::Bool

    z_primal::Tv
    z_dual::Tv

    x::Vector{Tv}
    Ax::Vector{Tv}
    y_lower::Vector{Tv}
    y_upper::Vector{Tv}
    s_lower::Vector{Tv}
    s_upper::Vector{Tv}

    Solution{Tv}(m, n) where{Tv} = new{Tv}(
        m, n, Sln_Unknown, Sln_Unknown, false, false,
        zero(Tv), zero(Tv),
        zeros(Tv, n), zeros(Tv, m),
        zeros(Tv, m), zeros(Tv, m),
        zeros(Tv, n), zeros(Tv, n)
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