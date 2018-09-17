using SparseArrays

function cholesky(A::AbstractSparseMatrix, d::AbstractVector) where{Ta<:Real}
    F = LinearAlgebra.cholesky(Symmetric(A*sparse(Diagonal(d))*A'))
    return F
end

function cholesky!(A::AbstractSparseMatrix, d::AbstractVector, F) where{Ta<:Real}
    # update Cholesky factor
    F = LinearAlgebra.cholesky!(F, Symmetric(A*sparse(Diagonal(d))*A'))
    return F
end