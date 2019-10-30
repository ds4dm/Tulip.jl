module TLPLinearAlgebra

using LinearAlgebra
BlasReal = LinearAlgebra.BlasReal

using SparseArrays

export construct_matrix

"""
    TLPLinearSolver{T<:Real}

Abstract container for linear solver used in solving the augmented system.
"""
abstract type TLPLinearSolver{Tv<:Real} end


include("sparse.jl")
include("dense.jl")

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
    addcolumn!(A, c)

Add column `c` to matrix `A`.
"""
addcolumn!(A::AbstractMatrix, c::AbstractVector) = (hcat(A, c), size(A, 2) + 1)

"""
    addrow!(A, r)

Add row `c` to matrix `A`.
"""
addrow!(A::AbstractMatrix, r::AbstractVector) = (vcat(A, r'), size(A, 1) + 1)

# include("unitBlockAngular.jl")


end  # module