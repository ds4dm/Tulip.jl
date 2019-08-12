module TLPLinearAlgebra

using LinearAlgebra
using SparseArrays

export factor_normaleq, factor_normaleq!, symbolic_cholesky, construct_matrix


"""
    construct_matrix

"""
function construct_matrix(
    ::Type{Matrix}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{Tv}
) where{Tv<:Real}

    A = zeros(Tv, m, n)
    for(i, j, v) in zip(aI, aJ, aV)
        A[i, j] = v
    end
    return A
end

construct_matrix(
    ::Type{SparseMatrixCSC}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{Tv}
) where{Tv<:Real} = sparse(aI, aJ, aV, m, n)



"""
    symbolic_cholesky
    
Compute Cholesky factorization of A*A'
"""
function symbolic_cholesky(A::AbstractMatrix{T}) where {T<:Real}

    F = factor_normaleq(A, ones(size(A, 2)))
    return F

end

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
    S = (A*Diagonal(d)*A')

    F = LinearAlgebra.cholesky(Symmetric(S), check=false)
    if !issuccess(F)
        
        # add regularization and try factor again.
        F = LinearAlgebra.cholesky!(F, Symmetric(S + 1e-6I), check=false)

        issuccess(F) || throw(PosDefException(2))
    end
    return F
end

function factor_normaleq!(
    A::AbstractMatrix{Tv},
    d::Vector{Tv},
    F::Factorization{Tv}
) where{Tv<:Real}
    # update Cholesky factor
    S = (A*Diagonal(d)*A')

    F = LinearAlgebra.cholesky!(F, Symmetric(S), check=false)
    if !issuccess(F)
        
        # add regularization and try factor again.
        F = LinearAlgebra.cholesky!(F, Symmetric(S + 1e-6I), check=false)

        # If factorization failed, throw error
        issuccess(F) || throw(PosDefException(2))
    end
    return F
end


function factor_normaleq!(A::Matrix{Tv}, d::Vector{Tv}, F_::Factorization{Tv}) where{Tv<:Real}
    S = (A*Diagonal(d)*A')
    
    F = LinearAlgebra.cholesky!(Symmetric(S), check=false)
    if !issuccess(F)
        
        # add regularization and try factor again.
        F = LinearAlgebra.cholesky!(Symmetric(S + 1e-6I), check=false)

        # If factorization failed, throw error
        issuccess(F) || throw(PosDefException(2))
    end
    return F
end


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

include("unitBlockAngular.jl")

end  # module