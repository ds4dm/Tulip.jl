function cholesky(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}) where{Ta<:Real, Td<:Real}
    F = cholfact(Symmetric(A*spdiagm(d)*A'))
    return F
end

function cholesky!(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}, F::Base.SparseArrays.CHOLMOD.Factor{Ta}) where{Ta<:Real, Td<:Real}
    # update Cholesky factor
    F = cholfact!(F, Symmetric(A*spdiagm(d)*A'))
    return F
end