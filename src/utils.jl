# Positive and negative part of a number
pos_part(x::T) where{T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where{T} = x >= zero(T) ? zero(T) : -x


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
