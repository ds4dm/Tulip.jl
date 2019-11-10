module TLPLinearAlgebra

using LinearAlgebra
BlasReal = LinearAlgebra.BlasReal

using SparseArrays

export construct_matrix

"""
    AbstractLinearSolver{Tv, Ta}

Abstract container for linear solver used in solving the augmented system.
"""
abstract type AbstractLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} end

# TODO: Traits

"""
    IndefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the indefinite augmented system.
"""
abstract type IndefLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} <: AbstractLinearSolver{Tv, Ta} end

"""
    PosDefLinearSolver{Tv, Ta}

Abstract container for linear solver working on the PSD normal equations system.
"""
abstract type PosDefLinearSolver{Tv<:Real, Ta<:AbstractMatrix{Tv}} <: AbstractLinearSolver{Tv, Ta} end


"""
    construct_matrix(Ta, m, n, aI, aJ, aV)

Construct matrix given matrix type `Ta`, size `m, n`, and data in COO format.
"""
function construct_matrix end

include("dense.jl")
include("sparse.jl")

# TODO: use parameter to choose between Indef/PosDef system
AbstractLinearSolver(A::Matrix{Tv}) where{Tv<:Real} = DenseLinearSolver(A)
AbstractLinearSolver(A::SparseMatrixCSC{Tv, Int64}) where{Tv<:BlasReal} = SparsePosDefLinearSolver(A)

end  # module