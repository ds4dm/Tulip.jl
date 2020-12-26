module TLPLinearAlgebra

using LinearAlgebra
using SparseArrays
export construct_matrix
export ZeroMatrix

import ..Tulip.Factory

"""
    construct_matrix(Ta, m, n, aI, aJ, aV)

Construct matrix given matrix type `Ta`, size `m, n`, and data in COO format.
"""
function construct_matrix end

function construct_matrix(
    ::Type{Matrix}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{T}
) where{T}
    A = zeros(T, m, n)
    # TODO: may be more efficient to first sort indices so that
    # A is accessed in column-major order.
    for(i, j, v) in zip(aI, aJ, aV)
        A[i, j] = v
    end
    return A
end

construct_matrix(
    ::Type{SparseMatrixCSC}, m::Int, n::Int,
    aI::Vector{Int}, aJ::Vector{Int}, aV::Vector{T}
) where{T} = sparse(aI, aJ, aV, m, n)

include("zero.jl")

end  # module
