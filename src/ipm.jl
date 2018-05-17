# import Tulip.Cholesky:
#     AbstractCholeskyFactor, SimpleDenseCholeskyFactor, SimpleSparseCholeskyFactor

"""
    TulipSolution contains a primal-dual point
"""
struct PrimalDualPoint{Tx<:Real, Ty<:Real, Ts<:Real, Tw<:Real, Tz<:Real}

    x::AbstractVector{Tx}  # primal variables
    y::AbstractVector{Ty}  # dual variables
    s::AbstractVector{Ts}  # dual variables wrt/ non-negativity of s

    w::AbstractVector{Tw}  # primal slack wrt/ upper bounds
    z::AbstractVector{Tz}  # dual variables wrt/ non-negativity of w
end


"""
    AbstractTulipModel
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
    n == size(c, 1) || throw(DimensionMismatch("Expected c to have size $(n) but got $(size(c, 1))"))
    m == size(b, 1) || throw(DimensionMismatch("Expected b to have size $(m) but got $(size(b, 1))"))
    uind[end] <= n  || throw(DimensionMismatch("Model has $(n) vars but upper bound given for var $(uind[end])"))

    s = PrimalDualPoint(similar(c), similar(b), similar(c), similar(uind), similar(uval))

    model = Model(m, n, A, b, c, uind, uval, s, :Built)
    return model

end


"""
    solve(model, tol, verbose)

Solve model m using an infeasible predictor-corrector Interior-Point algorithm.

# Arguments
- `model`: concrete model
- `tol`: numerical tolerance
- `verbose`: 0 means no output, 1 display logs at each iteration
"""
function solve!(
    model::Model;
    tol::Float64 = 10.0^-8,
    verbose::Int = 0
)

    N_ITER_MAX = 10  # maximum number of IP iterations
    niter = 0  # number of IP iterations

    # TODO: pre-optimization stuff

    F = symbolic_cholesky(model.A)  # Symbolic factorization
    θ = zeros(model.sol.x)

    # compute starting point
    compute_starting_point!(model, F)

    # TODO: check stopping criterion for possible early termination

    # IPM log
    println(" Itn      Primal Obj        Dual Obj  Prim Inf  Dual Inf  UBnd Inf\n")

    # main loop
    while niter < N_ITER_MAX
        
        # I. Form and factor Newton System
        compute_newton!(
            model.A,
            model.sol.x,
            model.sol.s,
            model.sol.w,
            model.sol.z,
            model.uind,
            θ,
            F
        )

        # II. Compute and take step
        compute_next_iterate!(model, F)

        # III. Book-keeping + display log
        niter += 1
        if verbose == 1
            print(@sprintf("%4d", niter))  # iteration count
            print(@sprintf("%+18.8e", 0.0))  # TODO: primal objective
            print(@sprintf("%+16.8e", 0.0))  # TODO: dual objective
            print(@sprintf("%10.2e", 0.0))  # TODO: primal infeas
            print(@sprintf("%9.2e", 0.0))  # TODO: dual infeas
            print(@sprintf("%9.2e", 0.0))  # TODO: upper bound infeas
            print("\n")
        end

        # check status
        if model.status == :Optimal
            println("Optimal solution found.")
            return model.status
        end

    end

    # 

    return model.status
    
end


"""
    compute_starting_point!
    Compute a starting point
"""
function compute_starting_point!(model::Model, F::Factorization)
    warn("TODO: starting point implementation")
    return model.sol
end

function compute_next_iterate!(model::Model, F::Factorization)
    warn("Implement computation of next iterate!")
    return model.sol
end

"""
    symbolic_cholesky
    Compute Cholesky factorization of A*A'
"""
function symbolic_cholesky(A::AbstractMatrix{T}) where {T<:Real}

    F = cholfact(Symmetric(A*A'))
    return F

end


"""
    compute_newton!
    Form and factorize the Newton system, using the normal equations.
"""
function compute_newton!(
    A::AbstractMatrix{Ta},
    x::AbstractVector{Tx},
    s::AbstractVector{Ts},
    w::AbstractVector{Tw},
    z::AbstractVector{Tz},
    uind::AbstractVector{Ti},
    θ::AbstractVector{T},
    F::Factorization{Ta}
    ) where {Ta<:Real, Tx<:Real, Ts<:Real, Tw<:Real, Tz<:Real, Ti<:Integer, T<:Real}

    # Compute Θ = (X^{-1} S + W^{-1} Z)^{-1}
    θ = x ./ s
    # update coefficients for upper bounds
    for (i, j) in enumerate(uind)
        θ[j] = 1.0 / (s[j] / x[j] + z[i] / w[i])
    end

    # Form the normal equations matrix and compute its factorization
    Cholesky.cholesky!(A, θ, F)

    return θ
end


"""
    solve_newton
    Solve Newton system with the given right-hand side.
    Overwrites the input d
"""
function solve_newton!(
    A::AbstractMatrix{Ta},
    θ::AbstractVector{T1},
    F::Factorization{Ta},
    Λ::PrimalDualPoint,
    d::PrimalDualPoint,
    uind::AbstractVector{Ti},
    ξ_b::AbstractVector{T2},
    ξ_c::AbstractVector{T3},
    ξ_u::AbstractVector{T4},
    ξ_xs::AbstractVector{T5},
    ξ_wz::AbstractVector{T6},
) where {Ta<:Real, T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real, T6<:Real, Ti<:Integer}

    ξ_tmp = ξ_c - (ξ_xs ./ Λ.x)
    ξ_tmp[uind] += (ξ_wz - (Λ.z .* ξ_u)) ./ Λ.w

    d.y = F \ (ξ_b + A * (θ .* ξ_tmp))

    d.x = θ .* (A' * d.y - ξ_tmp)

    d.z = (Λ.z .* (-ξ_u + d.x[uind]) + ξ_wz) ./ Λ.w
    d.s = (ξ_xs - Λ.s .* d.x) ./ Λ.x
    d.w = (ξ_wz - Λ.w .* d.z) ./ Λ.z

    return d
end