module TLPLinearAlgebra

using LinearAlgebra
using SparseArrays

export factor_normaleq

"""
    factor_normal_eqn(A, d)

Compute a Cholesky factorization of `A*D*A'`, where `D=Diag(d)`.

    factor_normal_eq!(A, d, F)
    
Compute a Cholesky factorization of `A*D*A'`, where `D=Diag(d)`, and overwrite
`F` in the process.
"""
function factor_normaleq(
    A::AbstractMatrix,
    d::AbstractVector
) where{Ta<:Real}
    F = LinearAlgebra.cholesky(Symmetric(A*Diagonal(d)*A'))
    return F
end

function factor_normaleq!(
    A::AbstractMatrix,
    d::AbstractVector,
    F
) where{Ta<:Real}
    # update Cholesky factor
    F = LinearAlgebra.cholesky!(F, Symmetric(A*Diagonal(d)*A'))
    return F
end

function consolidate!(A::AbstractMatrix) end

"""
    addcolumn!(A, c)

Add column `c` to matrix `A`.
"""
addcolumn!(A::AbstractMatrix, c::AbstractVector) = (hcat(A, c), size(A, 2) + 1)

"""
    addrow!(A, r)

Add row `c` to matrix `A`.
"""
addrow!(A::AbstractMatrix, r::AbstractVector) = (vcat(A, r'), size(A, 1) + 1)

include("denseBlockAngular.jl")

end  # module