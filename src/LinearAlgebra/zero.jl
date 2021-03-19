import Base:
    getindex, size
import LinearAlgebra.mul!

"""
    ZeroMatrix{T} <: AbstractMatrix{T}

A type for representing matrices whose coefficients are all zero.
"""
struct ZeroMatrix{T} <: AbstractMatrix{T}
    m::Int
    n::Int

    ZeroMatrix{T}(m, n) where{T} = new{T}(m, n)
end

size(A::ZeroMatrix) = (A.m, A.n)
getindex(::ZeroMatrix{T}, ::Integer, ::Integer) where{T} = zero(T)

mul!(x::AbstractVector{T}, ::ZeroMatrix{T}, y::AbstractVector{T}, α::Number, β::Number) where{T} = rmul!(x, β)
