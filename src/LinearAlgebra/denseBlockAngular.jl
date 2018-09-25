# Base Interface
import Base:
    size, getindex, copy

import LinearAlgebra:
    mul!, ldiv!

# TODO: make sure this implementation is stable

"""
    DenseBlockAngular{Tv<:Real}

# Attributes
-`m`: Number of linking constraints
-`n`: Total number of columns
-`R`: Number of blocks
-`colptr`: Index of first column of each block. The columns of
    block `r` are columns `colptr[k] ... colptr[k+1]-1`.
-`cols`: List of R blocks of columns
"""
mutable struct DenseBlockAngular{Tv<:Real} <: AbstractMatrix{Tv}
    m::Int
    n::Int
    R::Int

    colptr::Vector{Int}
    blocks::Vector{Matrix{Tv}}
    B0::Matrix{Tv}
    B::Matrix{Tv}

    DenseBlockAngular(
        m::Ti, n::Ti, R::Ti,
        colptr::Vector{Ti},
        blocks::Vector{Matrix{Tv}},
        B0::Matrix{Tv},
        B::Matrix{Tv}
    ) where {Tv<:Real, Ti<:Integer} = new{Tv}(m, n, R, colptr, blocks, B0, B)
end

"""
    DenseBlockAngular(blocks, C)

Construct a dense block-angular matrix, given the list of lower blocks, 
and columns `C`.
"""
function DenseBlockAngular(blocks, C::AbstractMatrix{T}) where T<:Real

    R = size(blocks, 1)
    if R == 0
        # No blocks, early return
        return DenseBlockAngular(
            size(C, 1), size(C, 2), 0,
            Vector{Int}(undef, 0),
            Vector{Matrix{T}}(undef, 0),
            C,
            C
        )
    end

    m = size(C, 1)
    colptr = Vector{Int}(undef, R+2)
    colptr[1] = 1

    # Dimension check
    for r in Base.OneTo(R)
        m == size(blocks[r], 1) || throw(DimensionMismatch(
            "Blocks 1 and $r have inconsistent dimensions."
        ))
        colptr[r+1] = colptr[r] + size(blocks[r], 2)
    end
    colptr[R+2] = colptr[R+1] + size(C, 2)

    blocks_ = deepcopy(blocks)
    C_ = deepcopy(C)

    return DenseBlockAngular(m, colptr[end]-1, R, colptr, blocks_, C_, hcat(blocks_..., C_))
end

"""
    DenseBlockAngular(blocks)

"""
function DenseBlockAngular(blocks::Vector{Matrix{T}}) where T<:Real 
    if length(blocks) == 0
        return DenseBlockAngular(
            0, 0, 0,
            Vector{Int}(),
            Vector{Matrix{T}}(),
            Matrix{T}(0, 0),
            Matrix{T}(undef, 0, 0)
        )
    else
        return DenseBlockAngular(blocks, zeros(T, size(blocks[1], 1), 0))
    end
end

function consolidate!(A::DenseBlockAngular)
    A.B = hcat(A.blocks..., A.B0)
    return nothing
end


size(M::DenseBlockAngular) = (M.m+M.R, M.n)
function getindex(M::DenseBlockAngular{Tv}, i::Integer, j::Integer) where Tv
    if !(1 <= i <= (M.m+M.R) && 1 <= j <= (M.n)); throw(BoundsError()); end
    
    if i > M.R
        # find index of block that contains column j
        blockidx = searchsortedlast(M.colptr, j)
        if blockidx == M.R+1
            return M.B0[i-M.R, j - M.colptr[blockidx]+1]
        else
            # return corresponding coefficient
            return M.blocks[blockidx][i-M.R, j - M.colptr[blockidx]+1]
        end
    else
        # find if column j belongs to block i
        return (M.colptr[i] <= j < M.colptr[i+1]) ? oneunit(Tv) : zero(Tv)
    end
    
    # never reached, for type stability
    return zero(Tv)
end

copy(M::DenseBlockAngular) = DenseBlockAngular(
    M.m, M.n, M.R,
    deepcopy(M.colptr),
    deepcopy(M.blocks),
    deepcopy(M.B0),
    deepcopy(M.B)
)

# Matrix-vector Multiplication
function mul!(
    y::AbstractVector{Tv},
    A::DenseBlockAngular{Tv},
    x::AbstractVector{Tv}
) where{Tv<:Real}

    m, n = size(A)
    n == length(x) || throw(DimensionMismatch(
        "A has dimensions $(size(A)) but x has dimension $(length(x))")
    )
    m == length(y) || throw(DimensionMismatch(
        "A has dimensions $(size(A)) but y has dimension $(length(y))")
    )

    slicedsum!(y, A.R, A.colptr, x)
    
    # Lower part of `y` is a matrix-vector product
    @views mul!(y[(A.R+1):end], A.B, x)
    return y
end

# Old version, uses list of blocks
# function mul!(
#     y::AbstractVector{Tv},
#     A::DenseBlockAngular{Tv},
#     x::AbstractVector{Tv}
# ) where{Tv<:Real}

#     m, n = size(A)
#     n == length(x) || throw(DimensionMismatch(
#         "A has dimensions $(size(A)) but x has dimension $(length(x))")
#     )
#     m == length(y) || throw(DimensionMismatch(
#         "A has dimensions $(size(A)) but y has dimension $(length(y))")
#     )
    
#     y_ = view(y, (A.R+1):(A.R+A.m))
#     y_ .= zero(Tv)

#     for r in 1:A.R
#         x_ = view(x, A.colptr[r]:(A.colptr[r+1]-1))
#         y[r] = sum(x_)
#         BLAS.gemv!('N', 1.0, A.blocks[r], x_, 1.0, y_)
#     end
#     @views BLAS.gemv!('N', 1.0, A.B0, x[A.colptr[A.R+1]:end], 1.0, y_)
    
#     return y
# end

"""
    mul!(x, At, y)

Compute Matrix-vector product `A'*y` and overwrites the result in `x`.
"""
function mul!(
    x::AbstractVector{Tv},
    At::LinearAlgebra.Transpose{Tv, DenseBlockAngular{Tv}},
    y::AbstractVector{Tv}
) where{Tv<:Real}
    
    A = At.parent

    m, n = size(A)
    n == length(x) || throw(DimensionMismatch(
        "A has dimensions $(size(A)) but x has dimension $(length(x))")
    )
    m == length(y) || throw(DimensionMismatch(
        "A has dimensions $(size(A)) but y has dimension $(length(y))")
    )

    @views mul!(x, transpose(A.B), y[(A.R+1):end])

    r = 0
    N = A.colptr[A.R+1]-1
    for i in Base.OneTo(N)
        @inbounds if i > (A.colptr[r+1]-1)
            r += 1
        end
        @inbounds x[i] += y[r]
    end
    
    return x
end

# Old version, uses the list of blocks
# """
#     mul!(x, At, y)

# Compute Matrix-vector product `A'*y` and overwrites the result in `x`.
# """
# function mul!(
#     x::AbstractVector{Tv},
#     At::LinearAlgebra.Transpose{Tv, DenseBlockAngular{Tv}},
#     y::AbstractVector{Tv}
# ) where{Tv<:Real}
    
#     A = At.parent

#     m, n = size(A)
#     n == length(x) || throw(DimensionMismatch(
#         "A has dimensions $(size(A)) but x has dimension $(length(x))")
#     )
#     m == length(y) || throw(DimensionMismatch(
#         "A has dimensions $(size(A)) but y has dimension $(length(y))")
#     )
    
#     y_ = y[(A.R+1):end]
    
#     @inbounds for r in 1:A.R
#         x_ = view(x, A.colptr[r]:(A.colptr[r+1]-1))
#         x_ .= y[r]
#         BLAS.gemv!('T', 1.0, A.blocks[r], y_, 1.0, x_)
#     end
#     @views BLAS.gemv!('T', 1.0, A.B0, y_, 0.0, x[A.colptr[A.R+1]:end])
    
#     return x
# end


"""
    FactorBlockAngular

# Attributes
-`m`: number of linking constraints
-`n`: total number of columns
-`R`: number of blocks
-`colptr`: index of first column of each block. The columns of
    block `r` are columns `colptr[k] ... colptr[k+1]-1`.
-`d`: `R`-dimensional vector containing the diagonal coefficients
    of the implicit Cholesky factor.
-`η`: An `m x R` matrix containing the lower blocks of the
    implicit Cholesky factor.
-`Fc`: Cholesky factor of the dense Schur complement.
"""
mutable struct FactorBlockAngular{Tv<:Real} <: Factorization{Tv}
    m::Int
    n::Int
    R::Int

    colptr::Vector{Int}
    
    d::Vector{Tv}
    η::Matrix{Tv}
    
    Fc::LinearAlgebra.Cholesky{Tv, Matrix{Tv}}

    FactorBlockAngular(
        m::Ti, n::Ti, R::Ti,
        colptr::Vector{Ti}, d::Vector{Tv}, η::Matrix{Tv},
        Fc
    ) where{Tv<:Real, Ti<:Integer} = new{Tv}(m, n, R, colptr, d, η, Fc)
end

size(F::FactorBlockAngular) = (F.m+F.R, F.n)


"""
    ldiv!(F, b)

In-place computation of `F \\ b`, overwriting `b` with the result. 
"""
function ldiv!(F::FactorBlockAngular{Tv}, b::AbstractVector{Tv}) where{Tv<:Real}
    
    # Dimension check
    (F.m + F.R) == length(b) || throw(DimensionMismatch(
        "F has dimensions $(size(F)) but b has dimensions $(length(b))"
    ))
    
    @views b[1:F.R] .*= F.d
    
    # right-hand side for global solve (b0)
    @views BLAS.gemv!('N', -1.0, F.η, b[1:F.R], 1.0, b[(F.R+1):end])
    
    # global solve
    @views ldiv!(F.Fc, b[(F.R+1):end])
    
    # local backward pass
    @views b[1:F.R] ./= F.d
    @views BLAS.gemv!('T', -1.0, F.η, b[(F.R+1):end], 1.0, b[1:F.R])
    @views b[1:F.R] .*= F.d
    
    return b
end

"""
    ldiv!(y, F, b)

Compute `F \\ b` and overwrite `y` with the result. `F`, `b` are not modified.
"""
function ldiv!(
    y::AbstractVector{Tv},
    F::FactorBlockAngular{Tv},
    b::AbstractVector{Tv}
) where{Tv<:Real}
    
    # Dimension check
    (F.m + F.R) == length(b) || throw(DimensionMismatch(
        "F has dimensions $(size(F)) but b has dimensions $(length(b))"
    ))
    (F.m + F.R) == length(y) || throw(DimensionMismatch(
        "F has dimensions $(size(F)) but y has dimensions $(length(y))"
    ))
    
    # over-write y with b
    y .= b

    ldiv!(F, y)
    return y

    # copy!(y, b)
    # @views y[1:F.R] .*= F.d
    
    # # right-hand side for global solve
    # @views BLAS.gemv!('N', -1.0, F.η, y[1:F.R], 0.0, y[(F.R+1):end])
    # @views BLAS.axpy!(1.0, b[(F.R+1):end], y[(F.R+1):end])
    
    # # global solve
    # @views ldiv!(F.Fc, y[(F.R+1):end])
    
    # # local backward pass
    # @views BLAS.gemv!('T', -1.0, F.η, y[(F.R+1):end], 0.0, y[1:F.R])
    # @views BLAS.axpy!(1.0, b[1:F.R], y[1:F.R])
    # @views y[1:F.R] .*= F.d
    
    # return y
end



"""
    CpBDBt!(C, B, d)

Compute C += B*D*B'
    where B' is the transpose of B, D = diag(d1, ..., dn),
    and B is mxn dense, C is mxm (dense), d is nx1 (dense).

C is modified in-place (no memory allocated),
only its upper-triangular part is modified.

# Arguments
-`C`: An `m`-by-`m` dense matrix, modified in place
-`B`: An `m`-by-`k` dense matrix
-`d`: A vector of length `k`
"""
function CpBDBt!(C::StridedMatrix, B::StridedMatrix, d::StridedVector)
    # TODO: allow user to specify which triangular part should be updated
    # TODO: add parameters for computing C = α C + β B*D*B'
    
    # Dimension checks
    (m, k) = size(B)
    (m, m) == size(C) || throw(DimensionMismatch(
        "C has dimensions $(size(C)) but B has dimensions $(size(B))"
    ))
    k == length(d) || throw(DimensionMismatch(
        "d has dimensions $(length(d)) but B has dimensions $(size(B))"
    ))

    # Linear scaling + BLAS call
    B_ = B * Diagonal(sqrt.(d))
    BLAS.syrk!('U', 'N', 1.0, B_, 1.0, C)
    
    return C
end

function factor_normaleq!(
    A::DenseBlockAngular{Tv},
    θ::AbstractVector{Tv},
    F::FactorBlockAngular{Tv}
    ) where{Tv<:Real}
    
    # Schur complement
    C = zeros(Tv, A.m, A.m)
    CpBDBt!(C, A.B, θ)

    # Diagonal of Cholesky factor
    slicedsum!(F.d, A.R, A.colptr, θ)
    F.d .\= oneunit(Tv)

    # Compute lower block of Cholesky factor
    for r in Base.OneTo(A.R)
        @views mul!(F.η[:, r], A.blocks[r], θ[A.colptr[r]:(A.colptr[r+1]-1)])
    end

    # low-rank update of C
    C_ = zeros(Tv, A.m, A.m)
    CpBDBt!(C_, F.η, F.d)
    C .-= C_

    # Cholesky factorization of (dense) Schur complement
    F.Fc = LinearAlgebra.cholesky(Symmetric(C), Val(false))
    
    return F
end

function factor_normaleq(
    A::DenseBlockAngular{Tv},
    θ::AbstractVector{Tv}
    ) where{Tv<:Real}
    
    # Schur complement
    C = zeros(Tv, A.m, A.m)
    d = zeros(Tv, A.R)
    η = zeros(Tv, A.m, A.R)
    
    # Compute schur-complement
    CpBDBt!(C, A.B, θ)

    # Compute diagonal part
    slicedsum!(d, A.R, A.colptr, θ)
    d .\= oneunit(Tv)

    # Compute lower block of Cholesky factor
    for r in Base.OneTo(A.R)
        @views mul!(η[:, r], A.blocks[r], θ[A.colptr[r]:(A.colptr[r+1]-1)])
    end

    # low-rank update of C
    C_ = zeros(Tv, A.m, A.m)
    CpBDBt!(C_, η, d)
    C .-= C_
    
    # Cholesky factorization of (dense) Schur complement
    S = Symmetric(C)
    Fc = LinearAlgebra.cholesky(S, Val(false))
    
    F = FactorBlockAngular(A.m, A.n, A.R, A.colptr, d,η, Fc)
end

"""
    addcolumn!(A, a)

Add column `a` to PrimalBlockAngular matrix `A`.
The index of the block is defined as the index of the first non-zero
element of `a[1:A.R]`.
"""
function addcolumn!(
    A::DenseBlockAngular{Tv},
    a::AbstractVector{Tv}
) where Tv<:Real

    # Dimension check
    (A.m + A.R) == length(a) || throw(DimensionMismatch(
        "A has $(A.m + A.R) rows but a has dimension $(length(a))"
    ))

    # Find index of block to which column pertains,
    # i.e. the smallest index 1 <= i <= A.R such that `a[i]` is non-zero.
    # If no such index is found, the column is added as a linking column.
    blockidx = A.R+1
    for r in Base.OneTo(A.R)
        if !iszero(a[r])
            blockidx = r
            break
        end
    end

    addcolumn!(A, a[(A.R+1):end], blockidx)
end

"""
    addcolumn!(A, a, i)

Add column `a` to block `i` of matrix `A`.
"""
function addcolumn!(
    A::DenseBlockAngular{Tv},
    a::AbstractVector{Tv},
    i::Int
) where Tv<:Real

    # Dimension check
    A.m == length(a) || throw(DimensionMismatch(""))

    if 1 <= i <= A.R
        # add column
        A.blocks[i] = hcat(A.blocks[i], a)
        k = A.colptr[i+1]

    elseif i == (A.R+1)
        A.B0 = hcat(A.B0, a)
        k = A.n + 1
    else
        throw(DimensionMismatch(
            "Attempting to add column to block $(i) while A has $(A.R) blocks"
        ))
    end
        
    # book-keeping
    A.n += 1
    for r in (i+1):(A.R+2)
        A.colptr[r] += 1
    end
    return A, k
end

function block_update!(
    r::Int,
    A::StridedMatrix{Tv},
    θ::StridedVector{Tv},
    η::StridedVector{Tv},
    d::StridedVector{Tv},
    C::StridedMatrix{Tv}
) where Tv<:Real
    
    d_ = oneunit(Tv) / sum(θ)
    BLAS.gemv!('N', 1.0, A, θ, 0.0, η)  # η = A * θ
    CpBDBt!(C, A, θ)    # C += 
    BLAS.syr!('U', -d_, η, C)   # C -= 

    d[r] = d_

    return nothing

end

"""
    slicedsum!(s, R, colptr, x)

Compute a sliced sum of `x` and overwrite `s` with the result
"""
function slicedsum!(s, R, colptr, x)
    length(s) >= R || throw(DimensionMismatch("s has wrong dimensions"))
    length(colptr) > R || throw(DimensionMismatch("colptr wrong dim"))
    @assert length(x) >= colptr[R+1]-1

    for r in Base.OneTo(R)
        @inbounds s[r] = zero(eltype(s))
        for i in UnitRange(colptr[r], colptr[r+1]-1)
            @inbounds s[r] += x[i]
        end
    end
    return s
end
