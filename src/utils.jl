using CodecBzip2
using CodecZlib

"""
    _open(f, fname)

Open a file with decompression stream as required.
"""
function _open(f::Function, fname::String)

    ext = Symbol(split(fname, ".")[end])

    if ext == :gz
        return Base.open(f, CodecZlib.GzipDecompressorStream, fname, "r")
    elseif ext == :bz2
        return Base.open(f, CodecBzip2.Bzip2DecompressorStream, fname, "r")
    else
        return Base.open(f, fname, "r")
    end
end

# Positive and negative part of a number
pos_part(x::T) where{T} = !signbit(x) ? x : zero(T)
neg_part(x::T) where{T} = !signbit(x) ? zero(T) : -x


@inline tones(Tv, n)  = fill!(Tv(undef, n),  one(eltype(Tv)))
@inline tzeros(Tv, n) = fill!(Tv(undef, n), zero(eltype(Tv)))

"""
    Factory{T}

Factory-like struct for passing options to lower-level components.
"""
struct Factory{T}
    T::Type{T}
    options::Base.Iterators.Pairs

    # Constructors
    Factory(::Type{T}; kwargs...) where{T} = new{T}(T, kwargs)
    Factory{T}(;kwargs...) where{T} = new{T}(T, kwargs)
end

instantiate(f::Factory{T}, args...; kwargs...) where{T} = T(args...; kwargs..., f.options...)
