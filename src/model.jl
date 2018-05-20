import Base: copy, +
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

    x::AbstractVector{T}  # primal variables
    w::AbstractVector{T}  # primal slack wrt/ upper bounds

    y::AbstractVector{T}  # dual variables
    s::AbstractVector{T}  # dual variables wrt/ non-negativity of s
    z::AbstractVector{T}  # dual variables wrt/ non-negativity of w
end

# extend Base operation on primal-dual points
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
"""
mutable struct Model{T1<:Real, T2<:Real, T3<:Real, T4<:Real, Ti<:Integer}

    # Problem Data
    nconstr::Integer  # number of constraints
    nvars::Integer  # number of variables
    
    A::AbstractMatrix{T1}  # constraint matrix
    b::AbstractVector{T2}  # right-hand side
    c::AbstractVector{T3}  # objective
    uind::AbstractVector{Ti}  # indices of variables with upper bounds
    uval::AbstractVector{T4}  # values of upper bounds on primal variables

    # current solution
    sol::PrimalDualPoint
    status::Symbol  # optimization status

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