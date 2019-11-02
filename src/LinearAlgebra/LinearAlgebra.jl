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


"""
    construct_matrix(Ta, m, n, aI, aJ, aV)

Construct matrix given matrix type `Ta`, size `m, n`, and data in COO format.
"""
function construct_matrix end

include("sparse.jl")
include("dense.jl")

end  # module