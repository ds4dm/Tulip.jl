function cholesky(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}) where{Ta<:Real, Td<:Real}
    F = cholfact(Symmetric(A*spdiagm(d)*A'))
    return F
end

function cholesky!(A::SparseMatrixCSC{Ta, Int64}, d::AbstractVector{Td}, F::Base.SparseArrays.CHOLMOD.Factor{Ta}) where{Ta<:Real, Td<:Real}
    # update Cholesky factor
    F = cholfact!(F, Symmetric(A*spdiagm(d)*A'))
    return F
end

"""
    addcolumn!

Add column to the problem and return the index of the new variable. The new
variable is assumed to be non-negative. `A`, `b`, `c`, `uind` and `uval` are
modified.

# Arguments:
- `A::AbstractMatrix{Tv, Ti}`: the constraint matrix
- `b::AbstractVector{Tv}`: right-hand side of linear constraints
- `c::AbstractVector{Tv}`: objective coefficients
- `uind::AbstractVector{Ti}`: indices of upper-bounded variables (sorted)
- `uval::AbstractVector{Tv}`: upper bounds' values
- `a::AbstractVector{Tv}`: the column to be added
- `objcoeff::Tv`: objective coefficient of the new variable
- `u::Tv`: upper bound of the new variable (non-negative)
"""
addcolumn!(A::SparseMatrixCSC{Tv, Ti},a::AbstractVector{Tv}) where{Tv<:Real, Ti<:Integer} = (hcat(A, a), size(A, 2) + 1)