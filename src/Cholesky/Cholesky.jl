module Cholesky

export DenseBlockAngular, FactorBlockAngular

"""
    cholesky(A, d)
    Form the matrix S = A*Diag(d)*A' and computes its Cholesky factorization.
    The factorization overwrites F.
"""
function cholesky(A::AbstractMatrix{Tv}, d::AbstractVector{Td}) where{Tv<:Real, Td<:Real}
    warn("cholesky! must be implemented for $(typeof(A))")
    F = cholfact(Symmetric(A*spdiagm(d)*A'))
    return F
end

function cholesky(A::Matrix{Ta}, d::AbstractVector{Td}) where{Ta<:Real, Td<:Real}
    F = cholfact(Symmetric(A*Diagonal(d)*A'))
    return F
end

function cholesky(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}) where{Ta<:Real, Td<:Real}
    F = cholfact(Symmetric(A*spdiagm(d)*A'))
    return F
end


"""
    cholesky!(A, d, F)
    Form the matrix S = A*Diag(d)*A' and computes its Cholesky factorization.
    The factorization overwrites F.
"""
function cholesky!(A::AbstractMatrix{Ta}, d::AbstractVector{Td}, F::Factorization{Ta}) where{Ta<:Real, Td<:Real}
    warn("cholesky! must be implemented for $(typeof(A))")
    return F
end

function cholesky!(A::Matrix{Ta}, d::AbstractVector{Td}, F::Base.LinAlg.Cholesky{Ta,Matrix{Ta}}) where{Ta<:Real, Td<:Real}
    F = cholfact!(F, A*Diagonal(d)*A')
    return F
end

function cholesky!(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}, F::Base.SparseArrays.CHOLMOD.Factor{Ta}) where{Ta<:Real, Td<:Real}
    # update CHolesky factor
    F = cholfact!(F, Symmetric(A*spdiagm(d)*A'))
    return F
end

"""
    solve!(y, F, b)
    Solves the system F*y = b and overwrites y.
"""
function solve!(y::AbstractVector{T1}, F::Factorization{T2}, b::AbstractVector{T3}) where {T1<:Real, T2<:Real, T3<:Real}
    warn("Implement solve! for concrete type implementations!")
    return y
end

# include("denseCholesky.jl")
# include("sparseCholesky.jl")
include("denseBlockAngular.jl")


end  # module