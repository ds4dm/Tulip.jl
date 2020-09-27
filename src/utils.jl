# Positive and negative part of a number
pos_part(x::T) where{T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where{T} = x >= zero(T) ? zero(T) : -x


@inline tones(Tv, n)  = fill!(Tv(undef, n),  one(eltype(Tv)))
@inline tzeros(Tv, n) = fill!(Tv(undef, n), zero(eltype(Tv)))