module Cholesky

# export AbstractCholeskyFactor, SimpleDenseCholeskyFactor, SimpleSparseCholeskyFactor

"""
    cholesky!(A, d, M)
    Form the matrix S = A*Diag(d)*A' and computes its Cholesky factorization.
    The factorization overwrites F.
"""
function cholesky!(A::AbstractMatrix{Ta}, d::AbstractVector{Td}, F::Factorization{Ta}) where{Ta<:Real, Td<:Real}
    warn("cholesky! must be implemented for $(typeof(A))")
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


end  # module