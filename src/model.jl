# extend Base operation on primal-dual points
import Base:
    copy, +

"""
    PrimalDualPoint

# Attributes
- `x::AbstractVector{T}`: vector of primal variables
- `w::AbstractVector{T}`: vector of primal slacks with respect to upper bounds
- `y::AbstractVector{T}`: vector of dual variables
- `s::AbstractVector{T}`: vector of dual variables (reduced costs of `x`)
- `z::AbstractVector{T}`: vector of dual variables (reduced costs of `w`)
"""
mutable struct PrimalDualPoint{T<:Real}
    x::AbstractVector{T}
    w::AbstractVector{T}

    y::AbstractVector{T}
    s::AbstractVector{T} 
    z::AbstractVector{T}
end

copy(p::PrimalDualPoint) = PrimalDualPoint(copy(p.x), copy(p.w), copy(p.y), copy(p.s), copy(p.z))

+(p1::PrimalDualPoint, p2::PrimalDualPoint) = PrimalDualPoint(
    p1.x + p2.x,
    p1.w + p2.w,
    p1.y + p2.y,
    p1.s + p2.s,
    p1.z + p2.z
)

"""
    Model
    Data structure for a model

# Attributes
-`nconstr::Integer`: Number of constraints
-`nvars::Integer`: Number of variables
-`A::AbstractMatrix{T1<:Real}`: Constraint matrix
-`b::AbstractVector{T2<:Real}`: Right-hand side of the equality constraints
-`c::AbstractVector{T3<:Real}`: Objective coefficient
-`uind::AbstractVector{Ti<:Integer}`: Indices of upper-bounded variables
-`uval::AbstractVector{T4<:Real}`: Upper bounds on the variables. Only finite
    upper bounds are stored.
-`sol::PrimalDualPoint`: Current solution to the problem (may be infeasible at
    the beginning of the optimization)
-`status::Symbol`: Optimization status
"""
mutable struct Model{T1<:Real, T2<:Real, T3<:Real, T4<:Real, Ti<:Integer}

    nconstr::Integer
    nvars::Integer
    
    A::AbstractMatrix{T1}
    b::AbstractVector{T2}
    c::AbstractVector{T3}
    uind::AbstractVector{Ti}
    uval::AbstractVector{T4}

    sol::PrimalDualPoint
    status::Symbol

end

function Model(
        A::AbstractMatrix{T1},
        b::AbstractVector{T2},
        c::AbstractVector{T3},
        uind::AbstractVector{Ti},
        uval::AbstractVector{T4}
    ) where{T1<:Real, T2<:Real, T3<:Real, T4<:Real, Ti<:Integer}
    
    (m, n) = size(A)
    p = size(uind, 1)
    n == size(c, 1) || throw(DimensionMismatch("Expected c to have size $(n) but got $(size(c, 1))"))
    m == size(b, 1) || throw(DimensionMismatch("Expected b to have size $(m) but got $(size(b, 1))"))
    if p > 0
        uind[end] <= n  || throw(DimensionMismatch("Model has $(n) vars but upper bound given for var $(uind[end])"))
    end
    s = PrimalDualPoint(ones(n), ones(p), zeros(m), ones(n), ones(p))

    model = Model(m, n, A, b, c, uind, uval, s, :Built)
    return model

end