# Right now this is just `BigFloat`, but in principle it could be expanded to a whitelist
# that would include other mutable types.
const SupportedMutableArithmetics = BigFloat

buffer_for_dot_product(::Type{V}) where {V <: AbstractVector{<:Real}} =
    buffer_for(LinearAlgebra.dot, V, V)

buffer_for_dot_product(::Type{F}) where {F <: Real} =
    buffer_for_dot_product(Vector{F})

buffered_dot_product_to!(
    buf::B,
    result::F,
    x::V,
    y::V,
) where {B <: Any, F <: SupportedMutableArithmetics, V <: AbstractVector{F}} =
    buffered_operate_to!(buf, result, LinearAlgebra.dot, x, y)

function buffered_dot_product!!(
    buf::B,
    x::V,
    y::V,
) where {B <: Any, F <: SupportedMutableArithmetics, V <: AbstractVector{F}}
    ret = zero(F)
    ret = buffered_dot_product_to!(buf, ret, x, y)
    return ret
end

buffered_dot_product!!(::Nothing, x::V, y::V) where {F <: Real, V <: AbstractVector{F}} =
    dot(x, y)

struct DotWeightedSumBuffer{F <: Real, DotBuffer <: Any}
    tmp::F
    dot::DotBuffer

    function DotWeightedSumBuffer{F}() where {F <: Real}
        dot_buffer = buffer_for_dot_product(F)
        return new{F, typeof(dot_buffer)}(zero(F), dot_buffer)
    end
end

struct DotWeightedSumBufferDummy
    dot::Nothing

    DotWeightedSumBufferDummy() = new(nothing)
end

buffer_for_dot_weighted_sum(::Type{F}) where {F <: SupportedMutableArithmetics} =
    DotWeightedSumBuffer{F}()

buffer_for_dot_weighted_sum(::Type{F}) where {F <: Real} =
    DotWeightedSumBufferDummy()

function buffered_dot_weighted_sum_to_inner!(
    buf::DotWeightedSumBuffer{F},
    sum::F,
    vecs::NTuple{n, NTuple{2, <:AbstractVector{F}}},
    weights::NTuple{n, <:Real},
) where {n, F <: SupportedMutableArithmetics}
    sum = zero!!(sum)

    for i in 1:n
        weight = weights[i]
        (x, y) = vecs[i]

        buffered_dot_product_to!(buf.dot, buf.tmp, x, y)
        mul!!(buf.tmp, weight)

        sum = add!!(sum, buf.tmp)
    end

    return sum
end

buffered_dot_weighted_sum_to!(
    buf::DotWeightedSumBuffer{F},
    sum::F,
    vecs::NTuple{n, NTuple{2, <:AbstractVector{F}}},
    weights::NTuple{n, Int}) where {n, F <: SupportedMutableArithmetics} =
    # It seems like the specialization
    #  *(x::BigFloat, c::Int8)
    # could be more efficient than
    #  *(x::BigFloat, c::Int)
    # MPFR has separate functions for those, and Julia uses them,
    # there must be a good (performance) reason for that.
    buffered_dot_weighted_sum_to_inner!(buf, sum, vecs, map(Int8, weights))

function buffered_dot_weighted_sum!!(
    buf::DotWeightedSumBuffer{F},
    vecs::NTuple{n, NTuple{2, <:AbstractVector{F}}},
    weights::NTuple{n, Int},
) where {n, F <: SupportedMutableArithmetics}
    ret = zero(F)
    ret = buffered_dot_weighted_sum_to!(buf, ret, vecs, weights)
    return ret
end

buffered_dot_weighted_sum!!(
    buf::DotWeightedSumBufferDummy,
    vecs::NTuple{n, NTuple{2, <:AbstractVector{F}}},
    weights::NTuple{n, Int}) where {n, F <: Real} =
    mapreduce((vec2, weight) -> weight*dot(vec2...), +, vecs, weights, init = zero(F))
