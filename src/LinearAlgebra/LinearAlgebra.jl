module TLPLinearAlgebra

using LinearAlgebra
BlasReal = LinearAlgebra.BlasReal

using SparseArrays

export construct_matrix


"""
    construct_matrix(Ta, m, n, aI, aJ, aV)

Construct matrix given matrix type `Ta`, size `m, n`, and data in COO format.
"""
function construct_matrix end

include("LinearSolvers/LinearSolvers.jl")

end  # module