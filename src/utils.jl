# Positive and negative part of a number
pos_part(x::T) where{T} = x >= zero(T) ? x : zero(T)
neg_part(x::T) where{T} = x >= zero(T) ? zero(T) : -x