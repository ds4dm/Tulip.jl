# Base Interface
import Base:
    size, getindex, copy

import LinearAlgebra:
    mul!, ldiv!

import SparseArrays:
    sparse

"""
    UnitBlockAngular{Tv}

Unit block angular matrix, i.e., primal block angular matrix where diagonal
    blocks are unit rows.
"""
mutable struct UnitBlockAngular{Tv<:Real} <: AbstractMatrix{Tv}

    m::Int  # number of linking constraints
    n::Int  # number of non-linking columns
    n0::Int  # Number of linking columns

    M::Int  # Total number of constraints
    N::Int  # Total number of columns
    R::Int  # Number of blocks

    B::Matrix{Tv}          # Lower blocks, concatenated
    B_::Matrix{Tv}         # Copy of `B`, not to allocate memory several times
    blockidx::Vector{Int}  # Block index of each column in `B`
    B0::Matrix{Tv}         # Linking block
    B0_::Matrix{Tv}        # Copy of the linking block, same as `B_`


    function UnitBlockAngular(B::Matrix{Tv}, R, blockidx::Vector{Int}, B0::AbstractMatrix{Tv}) where Tv<:Real

        # Sanity checks
        m, n = size(B)
        m0, n0 = size(B0)

        m == m0 || throw(DimensionMismatch(
            "B has $(m) rows, but B0 has $(m0)"
        ))

        n == length(blockidx) || throw(DimensionMismatch(
            "B has $n column, but blockidx has $(length(blockidx))"
        ))

        imin, imax = extrema(blockidx)
        (0 <= imin <= imax <= R) || throw(DimensionMismatch(
            "Specified $R blocks but indices are in range [$imin, $imax]"
        ))

        A = new{Tv}()

        # Dimensions
        A.m  = m
        A.n  = n
        A.n0 = n0
        A.R  = R
        A.M  = A.R + A.m
        A.N  = A.n0 + A.n
    
        # Copy data to avoid modifying argument afterwards
        A.B = copy(B)
        A.B_ = copy(B)
        A.blockidx = copy(blockidx)
        A.B0 = copy(B0)

        return A
    end

end

mutable struct FactorUnitBlockAngular{Tv<:Real} <: Factorization{Tv}
    m::Int
    n0::Int
    n::Int

    R::Int
    blockidx::Vector{Int}
    
    d::Vector{Tv}
    η::Matrix{Tv}

    Fc::LinearAlgebra.Cholesky{Tv, Matrix{Tv}}

    FactorUnitBlockAngular(
        m, n0, n, R, blockidx, d::Vector{Tv}, η::Matrix{Tv}, Fc
    ) where Tv<:Real = new{Tv}(m, n0, n, R, copy(blockidx), copy(d), copy(η), copy(Fc))

end

# Useful constructors
function UnitBlockAngular(blocks::Vector{Matrix{Tv}}, B0::AbstractMatrix{Tv}) where Tv<:Real
    R = length(blocks)
    (m, n0) = size(B0)

    if R == 0
        return UnitBlockAngular(zeros(Tv, m, 0), R, Int[], B0)
    end
    
    B = hcat(blocks...)
    n = size(B, 2)

    blockidx = ones(Int64, n)

    # Sanity checks
    k = 1
    for r in Base.OneTo(R)
        n_ = size(blocks[r], 2)
        blockidx[k:(k+n_-1)] .= r
        k += n_
    end

    return UnitBlockAngular(B, R, blockidx, B0)
end

function UnitBlockAngular(blocks::Vector{Matrix{Tv}}) where Tv<:Real
    if length(blocks) == 0
        return UnitBlockAngular(blocks, zeros(Tv, 0, 0))
    else
        return UnitBlockAngular(blocks, zeros(Tv, size(blocks[1], 1), 0))
    end
end

# Base matrix interface
size(A::UnitBlockAngular) = (A.M, A.N)

function getindex(A::UnitBlockAngular{Tv}, i::Integer, j::Integer) where Tv<:Real
    
    # Sanity check
    ((1 <= i <= A.M) && (1 <= j <= A.N)) || throw(BoundsError())
 
    if (i <= A.R) && (j <= A.n0)
        return zero(Tv)
    
    elseif (i <= A.R) && (j > A.n0)
        # Check if column j belongs to block i
        if A.blockidx[j-A.n0] == i
            return oneunit(Tv)
        else
            return zero(Tv)
        end
    
    elseif (i > A.R) && (j <= A.n0)
        # Element is in linking block
        return A.B0[i-A.R, j]
    
    elseif (i > A.R) && (j > A.n0)
        # Element is in regular block
        return A.B[i-A.R, j-A.n0]
    
    else
        error("Index ($i, $j) not found")
    end
        
end

copy(A::UnitBlockAngular) = UnitBlockAngular(A.B, A.R, A.blockidx, A.B0)

"""
    sparse(A)

Efficient convertion to sparse matrix
"""
function sparse(A::UnitBlockAngular{Tv}) where Tv<:Real
    
    A_ = hcat(
        vcat(spzeros(A.R, A.n0), A.B0),
        vcat(
            sparse(A.blockidx, collect(1:A.n), ones(A.n), A.R, A.n),
            A.B
        )
    )
    return A_
end


# Custom low-level functions
"""
    slicedsum_!(y, R, blockidx, x)

Compute `y_{j} = sum{x_{i} | blockidx[i] == j}`
"""
function slicedsum_!(y, R::Int, blockidx::Vector{Int}, x)

    # Sanity check
    length(blockidx) == length(x) || throw(DimensionMismatch(
        "blockidx and x have different dimensions."
    ))
    length(y) == R || throw(DimensionMismatch(
        "R=$R but y has dimension $(length(y))"
    ))
    y .= zero(eltype(y))
    
    @inbounds for i in Base.OneTo(length(x))
        j = blockidx[i]
        y[j] += x[i]
    end
    return y
end

"""
    unit_lift!(y, R, blockidx, x)

Compute `∀i, x[i] += y[blokidx[i]]`.
"""
function unit_lift!(y, R::Int, blockidx::Vector{Int}, x)

    # Sanity checks
    length(blockidx) == length(x) || throw(DimensionMismatch(
        "blockidx and x have different dimensions."
    ))
    length(y) == R || throw(DimensionMismatch(
        "R=$R but y has dimension $(length(y))"
    ))

    @inbounds for i in Base.OneTo(length(x))
        j = blockidx[i]
        x[i] += y[j]
    end
    return x
end

# Base linear algebra
"""
    mul!(y, A, x)

In-place computation of `y = A*x`
"""
function mul!(
    y::Vector,
    A::UnitBlockAngular{Tv},
    x::Vector
) where Tv<:Real
    slicedsum_!(view(y, 1:A.R), A.R, A.blockidx, view(x, (A.n0+1):A.N))
    
    @views BLAS.gemv!('N', 1.0, A.B, x[(1+A.n0):end], 0.0, y[(A.R+1):A.M])
    @views BLAS.gemv!('N', 1.0, A.B0, x[1:A.n0], 1.0, y[(A.R+1):A.M])
    
    return y
end

"""
    mul!(x, At, y)

Compute Matrix-vector product `A'*y` and overwrites the result in `x`.
"""
function mul!(
    x::Vector,
    At::Union{
        LinearAlgebra.Transpose{Tv, UnitBlockAngular{Tv}},
        LinearAlgebra.Adjoint{Tv, UnitBlockAngular{Tv}}
    },
    y::Vector
) where{Tv<:Real}
    A = At.parent

    m, n = size(A)
    n == length(x) || throw(DimensionMismatch(
        "A has size $(size(A)) but x has size $(length(x))")
    )
    m == length(y) || throw(DimensionMismatch(
        "A has size $(size(A)) but y has size $(length(y))")
    )

    # `x = B' * y0`
    @views mul!(x[(A.n0+1):end], transpose(A.B), y[(A.R+1):end])
    @views mul!(x[1:A.n0], transpose(A.B0), y[(A.R+1):end])

    # `∀i, x[i] += y[blockidx[i]]`
    @views unit_lift!(y[1:A.R], A.R, A.blockidx, x[(A.n0+1):end])

    return x
end

"""
    rank_update!(α, C, β, B, B_, d)

Compute `C = α * C + β * (B*D*B')` where `D=Diag(d)`
"""
function rank_update!(α::Float64, C, β::Float64, B, B_, d::StridedVector)

    # Sanity checks
    (m, k) = size(B)
    size(B) == size(B_) || throw(DimensionMismatch(
        "B and B_ have different sizes"
    ))
    (m, m) == size(C) || throw(DimensionMismatch(
        "C has size $(size(C)) but B has size $(size(B))"
    ))
    k == length(d) || throw(DimensionMismatch(
        "B has size $(size(B)) but d has size $(size(d))"
    ))

    # Copy data into B_
    copyto!(B_, B)

    # Diagonal scaling
    rmul!(B_, Diagonal(sqrt.(d)))

    # Rank-k update
    BLAS.syrk!('U', 'N', β, B_, α, C)

    return C

end

"""
    rank_update!(α, C, β, B, d)

Compute `C = α * C + β * (B*D*B')` where `D=Diag(d)`
"""
function rank_update!(α::Float64, C, β::Float64, B, d::StridedVector)

    # Sanity checks
    (m, k) = size(B)
    (m, m) == size(C) || throw(DimensionMismatch(
        "C has size $(size(C)) but B has size $(size(B))"
    ))
    k == length(d) || throw(DimensionMismatch(
        "B has size $(size(B)) but d has size $(size(d))"
    ))

    # Copy data
    B_ = copy(B)

    # Diagonal scaling
    rmul!(B_, Diagonal(sqrt.(d)))

    # Rank-k update
    BLAS.syrk!('U', 'N', β, B_, α, C)

    return C

end

"""
    lowerfactor!(L, R, blockidx, B, d)

Compute lower block of Cholesky factor.
"""
function lowerfactor!(
    L::Matrix{Tv},
    R::Int,
    blockidx::Vector{Int},
    B::Matrix{Tv},
    d::StridedVector
) where Tv<:Real

    m, n = size(B)

    n == length(d) || throw(DimensionMismatch(
        "B has size $(size(B)) but d has size $(size(d))"
    ))
    (m, R) == size(L) || throw(DimensionMismatch(
        "L has wrong size"
    ))

    L .= zero(Tv)

    @inbounds for i in Base.OneTo(n)
        d_ = d[i]
        k = blockidx[i]

        @inbounds for j in Base.OneTo(m)
            L[j, k] += d_ * B[j, i]
        end

    end

    return L        
end

"""
    factor_normaleq(A, θ)

Compute Cholesky factorization of `A*Θ*A'`, where `Θ = Diag(θ)`.
"""
function factor_normaleq(A::UnitBlockAngular{Tv}, θ::Vector{Tv}) where Tv<:Real

    C = zeros(Tv, A.m, A.m)
    d = zeros(Tv, A.R)
    η = zeros(Tv, A.m, A.R)

    θ0 = view(θ, 1:A.n0)
    θ_ = view(θ, (1+A.n0):A.N)

    # Compute Schur complement
    rank_update!(zero(Tv), C, oneunit(Tv), A.B, A.B_, θ_)
    rank_update!(oneunit(Tv), C, oneunit(Tv), A.B0, θ0)

    # Diagonal of Cholesky factor
    slicedsum_!(d, A.R, A.blockidx, θ_)
    d .\= oneunit(Tv)

    # Compute lower block of Cholesky factor
    lowerfactor!(η, A.R, A.blockidx, A.B, θ_)

    # low-rank update of C
    rank_update!(oneunit(Tv), C, -oneunit(Tv), η, d)

    # Cholesky factor of Schur complement
    S = Symmetric(C)
    Fc = cholesky(S, Val(false))

    return FactorUnitBlockAngular(A.m, A.n0, A.n, A.R, A.blockidx, d, η, Fc)
end


"""
    factor_normaleq!(F, A, θ)

In-place computation of Cholesky factorization of `A*Θ*A'`, where `Θ = Diag(θ)`.
"""
function factor_normaleq!(
    A::UnitBlockAngular{Tv},
    θ::Vector{Tv},
    F::FactorUnitBlockAngular{Tv}
) where Tv<:Real

    θ0 = view(θ, 1:A.n0)
    θ_ = view(θ, (1+A.n0):A.N)
    C = zeros(Tv, A.m, A.m)

    # Compute Schur complement
    rank_update!(zero(Tv), C, oneunit(Tv), A.B, A.B_, θ_)
    rank_update!(oneunit(Tv), C, oneunit(Tv), A.B0, θ0)

    # Diagonal of Cholesky factor
    slicedsum_!(F.d, A.R, A.blockidx, θ_)
    F.d .\= oneunit(Tv)

    # Compute lower block of Cholesky factor
    lowerfactor!(F.η, A.R, A.blockidx, A.B, θ_)

    # low-rank update of C
    rank_update!(oneunit(Tv), C, -oneunit(Tv), F.η, F.d)

    # Cholesky factor of Schur complement
    S = Symmetric(C)
    F.Fc = cholesky(S, Val(false))

    return F
end

"""
    ldiv!(F, b)

In-place computation of `F \\ b`, overwriting `b` with the result.
"""
function ldiv!(
    F::FactorUnitBlockAngular{Tv},
    b::Vector
) where Tv<:Real

    (F.m+F.R) == length(b) || throw(DimensionMismatch(
        "F has size $(F.m+F.R) but b has size $(size(b))"
    ))

    @views b[1:F.R] .*= F.d

    @views BLAS.gemv!('N', -oneunit(Tv), F.η, b[1:F.R], oneunit(Tv), b[(F.R+1):end])
    
    # global solve
    @views ldiv!(F.Fc, b[(F.R+1):end])

    # local backward pass
    @views b[1:F.R] ./= F.d
    @views BLAS.gemv!('T', -oneunit(Tv), F.η, b[(F.R+1):end], oneunit(Tv), b[1:F.R])
    @views b[1:F.R] .*= F.d

    return b
end

"""
    ldiv!(y, F, b)

Compute `F \\ b` and store the result in `y`.
"""
function ldiv!(
    y::Vector{Tv},
    F::FactorUnitBlockAngular{Tv},
    b::Vector{Tv}
) where Tv<:Real

    # Dimension check
    (F.m + F.R) == length(b) || throw(DimensionMismatch(
        "F has dimensions $(F.m + F.R) but b has dimensions $(length(b))"
    ))
    (F.m + F.R) == length(y) || throw(DimensionMismatch(
        "F has dimensions $(F.m + F.R) but y has dimensions $(length(y))"
    ))
    y .= b

    ldiv!(F, y)
    return y
end