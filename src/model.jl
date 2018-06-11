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

# extend Base operation on primal-dual points
import Base:
    copy, +

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
-`n_var::Int`: Number of variables
-`n_con::Int`: Number of constraints
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
mutable struct Model
    #=======================================================
        Optimization environment
    =======================================================#

    env::TulipEnv


    #=======================================================
        Problem data
    =======================================================#

    n_var::Int              # Number of variables (total)
    n_var_ub::Int           # Number of upper-bounded variables
    n_con::Int              # Number of constraints

    A::AbstractMatrix       # Constraints matrix
    b::AbstractVector       # Right-hand side of equality constraints
    c::AbstractVector       # Objective vector
    uind::AbstractVector{Int}  # Indices of upper-bounded variables
    uval::AbstractVector    # Values of upper bounds


    #=======================================================
        Book-keeping
    =======================================================#

    sol::PrimalDualPoint    # Current primal-dual iterate
    status::Symbol          # Optimization status


    #=======================================================
        Model constructor
    =======================================================#  
    
    function Model(env, A, b, c, uind, uval, sol)
        m = new()

        m.env = env

        # Dimension check
        n_con = size(A, 1)
        n_var = size(A, 2)
        n_con == size(b, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but b has size $(size(b))"
        ))
        n_var == size(c, 1) || throw(DimensionMismatch(
            "A has size $(size(A)) but c has size $(size(b))"
        ))
        n_var_ub = size(uind, 1)
        n_var_ub == size(uval, 1) || throw(DimensionMismatch(
            "uind has size $(size(uval)) but uval has size $(size(uval))"
        ))
        n_var_ub <= n_var || throw(DimensionMismatch(
            "Too many upper bounds were specified"
        ))
        if n_var_ub > 0
            uind[end] <= n_var || throw(DimensionMismatch(
                "Got upper bound for var $(uind[end])>$(n_var)"
            ))
        end

        m.n_var = n_var
        m.n_var_ub = n_var_ub
        m.n_con = n_con
        m.A = A
        m.b = b
        m.c = c
        m.uind = uind
        m.uval = uval

        m.sol = sol
        m.status = :Built

        return m
    end
end

"""   
    Model(A, b, c, uind, uval)

Constructs a model with upper bounds on the specified variables

    Model(A, b, c)

Constructs a model with no upper bounds

    Model()

Empty model
"""
function Model(
    A::AbstractMatrix{T1},
    b::AbstractVector{T2},
    c::AbstractVector{T3},
    uind::AbstractVector{Ti},
    uval::AbstractVector{T4}
    ) where{T1<:Real, T2<:Real, T3<:Real, T4<:Real, Ti<:Integer}
    
    env = TulipEnv()

    (m, n) = size(A)
    p = size(uind, 1)

    sol0 = PrimalDualPoint(ones(n), ones(p), zeros(m), ones(n), ones(p))

    model = Model(
        env,
        A, b, c, uind, uval, sol0
    )
    return model

end

Model() = Model(spzeros(0, 0), Vector{Float64}(0,), Vector{Float64}(0,))

Model(A, b, c) = Model(A, b, c, Vector{Int}(0,), Vector{Float64}(0,))



