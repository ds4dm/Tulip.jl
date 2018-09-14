module LinearAlgebra

export DenseBlockAngular, FactorBlockAngular

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
include("sparseLinearAlgebra.jl")

end  # module