# Base Interface
import Base:
    size, getindex, copy
import Base.LinAlg:
    A_mul_B!, At_mul_B, A_ldiv_B!

# TODO: make sure this implementation is stable

"""
    DenseBlockAngular{Tv<:Real, Ti<:Integer}

# Attributes
-`m::Ti`: Number of linking constraints
-`n::Ti`: Total number of columns
-`R::Ti`: Number of blocks
-`colptr::Vector{Ti}`: Index of first column of each block. The columns of
    block `r` are columns `colptr[k] ... colptr[k+1]-1`.
-`cols::Vector{Matrix{Tv}}`: List of R blocks of columns
"""
mutable struct DenseBlockAngular{Tv<:Real} <: AbstractMatrix{Tv}
    m::Int
    n::Int
    R::Int

    colptr::Vector{Int}
    cols::Vector{Matrix{Tv}}
    colslink::Matrix{Tv}

    DenseBlockAngular(m::Ti, n::Ti, R::Ti, colptr::Vector{Ti}, cols::Vector{Matrix{Tv}}, colslink::Matrix{Tv}
    ) where {Tv<:Real, Ti<:Integer} = new{Tv}(m, n, R, colptr, cols, colslink)
end


function DenseBlockAngular(col_list::Vector{Matrix{T}}, colslink::Matrix{T}) where T<:Real

    R = size(col_list, 1)
    if R == 0
        return DenseBlockAngular(0, 0, 0, Vector{Int}(0), Matrix{T}(0, 0))
    end
    m = size(col_list[1], 1)

    colptr = Vector{Int}(R+2)
    colptr[1] = 1

    # Dimension check
    for r in 1:R
        m == size(col_list[r], 1) || throw(DimensionMismatch("Block 1 has dimensions $(size(col_list[1])) but block $(r) has dimensions $(size(col_list[r]))"))
        colptr[r+1] = colptr[r] + size(col_list[r], 2)
    end
    m == size(colslink, 1) || throw(DimensionMismatch("Block 1 has dimensions $(size(col_list[1])) but B has dimensions $(size(colslink))"))
    colptr[R+2] = colptr[R+1] + size(colslink, 2)

    cols = deepcopy(col_list)
    B = deepcopy(colslink)

    return DenseBlockAngular(m, colptr[end]-1, R, colptr, cols, B)
end
DenseBlockAngular(col_list::Vector{Matrix{T}}) where T<:Real = DenseBlockAngular(col_list, zeros(T, size(col_list[1], 1), 0))


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
            return M.cols[blockidx][i-M.R, j - M.colptr[blockidx]+1]
        end
    else
        # find if column j belongs to block i
        return (M.colptr[i] <= j < M.colptr[i+1]) ? oneunit(Tv) : zero(Tv)
    end
    
    # never reached, for type stability
    return zero(Tv)
end

copy(M::DenseBlockAngular) = DenseBlockAngular(M.m, M.n, M.R, deepcopy(M.colptr), deepcopy(M.cols), deepcopy(M.colslink))

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
        Base.BLAS.gemv!('N', 1.0, A.cols[r], x_, 1.0, y_)
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
        Base.BLAS.gemv!('T', 1.0, A.cols[r], y_, 1.0, x_)
    end
    @views Base.BLAS.gemv!('T', 1.0, A.colslink, y_, 1.0, x[A.colptr[A.R+1]:end])
    
    return x
end


"""
    FactorBlockAngular

# Attributes
-`m::Integer`: number of linking constraints
-`n::Integer`: total number of columns
-`R::Integer`: number of blocks
-`colptr::Vector{Int}`: index of first column of each block. The columns of
    block `r` are columns `colptr[k] ... colptr[k+1]-1`.
-`d::Vector{T}`: `R`-dimensional vector containing the diagonal coefficients
    of the implicit Cholesky factor.
-`η::Matrix{T}`: An `m x R` matrix containing the lower blocks of the
    implicit Cholesky factor.
-`Fc::Base.LinAlg.Cholesky{T, Matrix{T}}`: Cholesky factor of the dense
    Schur complement.
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
    where B' is the transpose of B, and D = diag(d1, ..., dn)
    and B is mxn dense, C is mxm (dense), d is nx1 (dense).

C is modified in-place (no memory allocated),
only its upper-triangular part is modified.

# Arguments
-`C::StridedMatrix{T}`: An `m`-by-`m` dense matrix, modified in place
-`B::StridedMatrix{T}`: An `m`-by-`k` dense matrix
-`d::StridedVector{T}`: A vector of length `k`
"""
function CpBDBt!(C::StridedMatrix{T}, B::StridedMatrix{T}, d::StridedVector{T}) where {T<:Real}
    # TODO: allow user to specify which triangular part should be update
    # TODO: add parameter for computing C = α C + B*D*B'
    
    # Dimension checks
    (m, k) = size(B)
    (m, m) == size(C) || throw(DimensionMismatch("C has dimensions $(size(C)) but B has dimensions $(size(B))"))
    k == size(d, 1) || throw(DimensionMismatch("d has dimensions $(size(d)) but B has dimensions $(size(B))"))
    
    # compute symmetric rank-k update
    # only the upper-triangular part of C is modified
    # temp = zero(T)
    # for j=1:m
    #     for l = 1:k
    #         @inbounds temp = d[l] * B[j, l]
    #         if temp == zero(T)
    #             # skip loop if temp is zero
    #             continue 
    #         end
    #         for i=1:j
    #             @inbounds C[i, j] = C[i, j] + temp * B[i, l]
    #         end
    #     end
    # end

    # # Linear scaling + BLAS call is more efficient
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
        # A_ = A.cols[:, A.colptr[r]:(A.colptr[r+1]-1)]
        η_ = view(F.η, :, r)
        
        # diagonal elements
        F.d[r] = oneunit(Tv) / sum(θ_)
        
        # lower factors
        Base.BLAS.gemv!('N', 1.0, A.cols[r], θ_, 0.0, η_)

        # Schur complement
        CpBDBt!(C, A.cols[r], θ_)
        Base.BLAS.syr!('U', -F.d[r], η_, C)
    end
    @views CpBDBt!(C, A.colslink, θ[A.colptr[(A.R+1)]:end])

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
        # A_ = A.cols[:, A.colptr[r]:(A.colptr[r+1]-1)]
        η_ = view(η, :, r)  # use view because η is modified in-place
        
        # diagonal elements
        d[r] = oneunit(Tv) / sum(θ_)
        
        # lower factors
        Base.BLAS.gemv!('N', 1.0, A.cols[r], θ_, 0.0, η_)
        
        # Schur complement
        CpBDBt!(C, A.cols[r], θ_)
        Base.BLAS.syr!('U', -d[r], η_, C)
    end
    @views CpBDBt!(C, A.colslink, θ[A.colptr[(A.R+1)]:end])
    
    # compute Cholesky factor of C
    Fc = cholfact(Symmetric(C))
    
    F = FactorBlockAngular(A.m, A.n, A.R, A.colptr, d,η, Fc)
end

function addcolumn!(A::DenseBlockAngular{Tv}, a::AbstractVector{Tv}, blockidx::Int) where Tv<:Real

    # Dimension check
    A.m == size(a, 1) || throw(DimensionMismatch(""))

    # add column
    A.cols[blockidx] = hcat(A.cols[blockidx], a)
    k = A.colptr[blockidx+1]

    # book-keeping
    A.n += 1
    for r in (blockidx+1):(A.R+2)
        A.colptr[r] += 1
    end

    return A, k
end

function addcolumn!(A::DenseBlockAngular{Tv}, a::AbstractVector{Tv}) where Tv<:Real

    # find index of block to which column pertains
    blockidx = 1
    for r in 1:A.R
        if a[r] > 0
            blockidx = r
            break
        end
    end

    addcolumn!(A, a[(A.R+1):end], blockidx)
end