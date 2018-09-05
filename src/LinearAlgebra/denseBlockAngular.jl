# Base Interface
import Base:
    size, getindex, copy
import Base.LinAlg:
    A_mul_B!, At_mul_B, A_ldiv_B!

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
    colslink::Matrix{Tv}

    DenseBlockAngular(m::Ti, n::Ti, R::Ti, colptr::Vector{Ti}, blocks::Vector{Matrix{Tv}}, colslink::Matrix{Tv}
    ) where {Tv<:Real, Ti<:Integer} = new{Tv}(m, n, R, colptr, blocks, colslink)
end

"""
    DenseBlockAngular(blocks, C)

Construct a dense block-angular matrix from a list of blocks and linking columns `C`.
"""
function DenseBlockAngular(blocks, C::AbstractMatrix{T}) where T<:Real

    R = size(blocks, 1)
    if R == 0
        # No blocks, early return
        return DenseBlockAngular(size(C, 1), size(C, 2), 0, Vector{Int}(0), Vector{Matrix{T}}(), C)
    end

    m = size(C, 1)
    colptr = Vector{Int}(R+2)
    colptr[1] = 1

    # Dimension check
    for r in Base.OneTo(R)
        m == size(blocks[r], 1) || throw(DimensionMismatch("Block 1 has dimensions $(size(col_list[1])) but block $(r) has dimensions $(size(col_list[r]))"))
        colptr[r+1] = colptr[r] + size(blocks[r], 2)
    end
    colptr[R+2] = colptr[R+1] + size(C, 2)

    cols = deepcopy(blocks)
    B = deepcopy(C)

    return DenseBlockAngular(m, colptr[end]-1, R, colptr, cols, B)
end

"""
    DenseBlockAngular(blocks)

"""
function DenseBlockAngular(blocks::Vector{Matrix{T}}) where T<:Real 
    if length(blocks) == 0
        return DenseBlockAngular(0, 0, 0, Vector{Int}(0), Vector{Matrix{T}}(), Matrix{T}(0, 0))
    else
        return DenseBlockAngular(blocks, zeros(T, size(blocks[1], 1), 0))
    end
end

size(M::DenseBlockAngular) = (M.m+M.R, M.n)
function getindex(M::DenseBlockAngular{Tv}, i::Integer, j::Integer) where {Tv<:Real}
    if !(1 <= i <= (M.m+M.R) && 1 <= j <= (M.n)); throw(BoundsError()); end
    
    if i > M.R
        # find index of block that contains column j
        blockidx = searchsortedlast(M.colptr, j)
        if blockidx == M.R+1
            return M.colslink[i-M.R, j - M.colptr[blockidx]+1]
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

copy(M::DenseBlockAngular) = DenseBlockAngular(M.m, M.n, M.R, deepcopy(M.colptr), deepcopy(M.blocks), deepcopy(M.colslink))

# Matrix-vector Multiplication
function A_mul_B!(y::AbstractVector{Tv}, A::DenseBlockAngular{Tv}, x::AbstractVector{Tv}) where{Tv<:Real}
    
    m, n = size(A)
    n == size(x, 1) || throw(DimensionMismatch("A has dimensions $(size(A)) but x has dimension $(size(x))"))
    m == size(y, 1) || throw(DimensionMismatch("A has dimensions $(size(A)) but y has dimension $(size(y))"))
    
    y_ = view(y, (A.R+1):(A.R+A.m))
    y_ .= zero(Tv)

    for r in 1:A.R
        x_ = view(x, A.colptr[r]:(A.colptr[r+1]-1))
        y[r] = sum(x_)
        Base.BLAS.gemv!('N', 1.0, A.blocks[r], x_, 1.0, y_)
    end
    @views Base.BLAS.gemv!('N', 1.0, A.colslink, x[A.colptr[A.R+1]:end], 1.0, y_)
    
    return y
end

function At_mul_B(A::DenseBlockAngular{Tv}, y::AbstractVector{Tv}) where{Tv<:Real}
    
    m, n = size(A)
    m == size(y, 1) || throw(DimensionMismatch("A has dimension $(size(A)) but y has dimension $(size(y))"))
    
    x = zeros(n)
    y_ = y[(A.R+1):end]
    
    @inbounds for r in 1:A.R
        x_ = view(x, A.colptr[r]:(A.colptr[r+1]-1))
        x_ .= y[r]
        Base.BLAS.gemv!('T', 1.0, A.blocks[r], y_, 1.0, x_)
    end
    @views Base.BLAS.gemv!('T', 1.0, A.colslink, y_, 1.0, x[A.colptr[A.R+1]:end])
    
    return x
end


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
    
    Fc::Base.LinAlg.Cholesky{Tv, Matrix{Tv}}

    FactorBlockAngular(
        m::Ti, n::Ti, R::Ti,
        colptr::Vector{Ti}, d::Vector{Tv}, η::Matrix{Tv},
        Fc::Base.LinAlg.Cholesky{Tv, Matrix{Tv}}
    ) where{Tv<:Real, Ti<:Integer} = new{Tv}(m, n, R, colptr, d, η, Fc)
end

size(F::FactorBlockAngular) = (F.m+F.R, F.n)


"""
    A_ldiv_B!

"""
function A_ldiv_B!(F::FactorBlockAngular{Tv}, b::AbstractVector{Tv}) where{Tv<:Real}
    
    # Dimension check
    (F.m + F.R) == size(b, 1) || throw(DimensionMismatch("F has dimensions $(size(F)) but b has dimensions $(size(b))"))
    
    @views b[1:F.R] .*= F.d
    
    # right-hand side for global solve (b0)
    @views Base.BLAS.gemv!('N', -1.0, F.η, b[1:F.R], 1.0, b[(F.R+1):end])
    
    # global solve
    @views Base.LinAlg.A_ldiv_B!(F.Fc, b[(F.R+1):end])
    
    # local backward pass
    @views b[1:F.R] ./= F.d
    @views Base.BLAS.gemv!('T', -1.0, F.η, b[(F.R+1):end], 1.0, b[1:F.R])
    @views b[1:F.R] .*= F.d
    
    return b
end

function A_ldiv_B!(y::AbstractVector{Tv}, F::FactorBlockAngular{Tv}, b::AbstractVector{Tv}) where{Tv<:Real}
    
    # Dimension check
    (F.m + F.R) == size(b, 1) || throw(DimensionMismatch("F has dimensions $(size(F)) but b has dimensions $(size(b))"))
    (F.m + F.R) == size(y, 1) || throw(DimensionMismatch("F has dimensions $(size(F)) but y has dimensions $(size(y))"))
    
    # over-write y with b
    copy!(y, b)
    @views y[1:F.R] .*= F.d
    
    # right-hand side for global solve
    @views Base.BLAS.gemv!('N', -1.0, F.η, y[1:F.R], 0.0, y[(F.R+1):end])
    @views Base.BLAS.axpy!(1.0, b[(F.R+1):end], y[(F.R+1):end])
    
    # global solve
    @views Base.LinAlg.A_ldiv_B!(F.Fc, y[(F.R+1):end])
    
    # local backward pass
    @views Base.BLAS.gemv!('T', -1.0, F.η, y[(F.R+1):end], 0.0, y[1:F.R])
    @views Base.BLAS.axpy!(1.0, b[1:F.R], y[1:F.R])
    @views y[1:F.R] .*= F.d
    
    return y
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
function CpBDBt!(C::StridedMatrix{T}, B::StridedMatrix{T}, d::StridedVector{T}) where {T<:Real}
    # TODO: allow user to specify which triangular part should be updated
    # TODO: add parameter for computing C = α C + B*D*B'
    
    # Dimension checks
    (m, k) = size(B)
    (m, m) == size(C) || throw(DimensionMismatch("C has dimensions $(size(C)) but B has dimensions $(size(B))"))
    k == size(d, 1) || throw(DimensionMismatch("d has dimensions $(size(d)) but B has dimensions $(size(B))"))

    # Linear scaling + BLAS call
    B_ = B * Diagonal(sqrt.(d))
    Base.BLAS.syrk!('U', 'N', 1.0, B_, 1.0, C)
    
    return C
end

# TODO: blocked version
# Compare to BLAS `syrk!`, there is a performance factor of ~2-3 to be gained 


function cholesky!(
    A::DenseBlockAngular{Tv},
    θ::AbstractVector{Tv},
    F::FactorBlockAngular{Tv}
    ) where{Tv<:Real}
    
    # Schur complement
    C = zeros(Tv, A.m, A.m)
    
    for r in 1:A.R
        # copying ararys uses more memory, but is faster
        θ_ = view(θ, A.colptr[r]:(A.colptr[r+1]-1))
        # A_ = A.blocks[:, A.colptr[r]:(A.colptr[r+1]-1)]
        η_ = view(F.η, :, r)
        
        block_update!(r, A.blocks[r], θ_, η_, F.d, C)
    end

    # Linking columns
    @views CpBDBt!(C, A.colslink, θ[A.colptr[(A.R+1)]:end])

    # Cholesky factorization of (dense) Schur complement
    F.Fc = cholfact(Symmetric(C))
    
    return F
end

function cholesky(
    A::DenseBlockAngular{Tv},
    θ::AbstractVector{Tv}
    ) where{Tv<:Real}
    
    # Schur complement
    C = zeros(Tv, A.m, A.m)
    d = zeros(Tv, A.R)
    η = zeros(Tv, A.m, A.R)
    
    for r in 1:A.R
        θ_ = view(θ, A.colptr[r]:(A.colptr[r+1]-1))
        η_ = view(η, :, r)  # use view because η is modified in-place
        
        block_update!(r, A.blocks[r], θ_, η_, d, C)
    end

    # Linking columns
    @views CpBDBt!(C, A.colslink, θ[A.colptr[(A.R+1)]:end])
    
    # Cholesky factorization of (dense) Schur complement
    Fc = cholfact(Symmetric(C))
    
    F = FactorBlockAngular(A.m, A.n, A.R, A.colptr, d,η, Fc)
end

"""
    addcolumn!(A, a)

"""
function addcolumn!(A::DenseBlockAngular{Tv}, a::AbstractVector{Tv}) where Tv<:Real

    # Dimension check
    (A.m + A.R) == length(a) || throw(DimensionMismatch("A has $(A.m + A.R) rows but a has dimension $(length(a))"))

    # Find index of block to which column pertains, i.e. the smallest index 1 <= i <= A.R such that `a[i]` is non-zero.
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

Add column `a` to block `i`. If 
"""
function addcolumn!(A::DenseBlockAngular{Tv}, a::AbstractVector{Tv}, i::Int) where Tv<:Real

    # Dimension check
    A.m == size(a, 1) || throw(DimensionMismatch(""))

    if 1 <= i <= A.R
        # add column
        A.blocks[i] = hcat(A.blocks[i], a)
        k = A.colptr[i+1]

    elseif i == (A.R+1)
        A.colslink = hcat(A.colslink, a)
        k = A.n + 1
    else
        throw(DimensionMismatch("Attempting to add to block $(i) while A has $(A.R) blocks"))
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
    Base.BLAS.gemv!('N', 1.0, A, θ, 0.0, η)
    CpBDBt!(C, A, θ)
    Base.BLAS.syr!('U', -d_, η, C)

    d[r] = d_

    return nothing

end