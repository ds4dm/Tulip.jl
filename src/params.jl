"""
    AbstractParam{T}

Abstract representation of a parameter with valuetype `T`.
"""
abstract type AbstractParam{T} end

"""
    get_param_name(::AbstractParam)

Return name of parameter.
"""
function get_param_name end

"""
    get_param_type(::AbstractParam)

Get parameter's value type.
"""
function get_param_type end

"""
    get_param_value(::AbstractParam)

Return current value of parameter.
"""
function get_param_value end

"""
    set_param_default!(::AbstractParam{T})

Set parameter to its default value.
"""
function set_param_default! end

"""
    set_param_value!(::AbstractParam{T}, v::T)

Check if `v` is an admissible value for parameter and, if so, change parameter's
    value to `v`.
"""
function set_param_value! end

"""
    test_param_value(::AbstractParam{T}, v::T)

Check whether value `v` is admissible for given parameter.
"""
function test_param_value end

"""
    RealParam{T<:Real}

Container for numerical (real-valued) parameters.
"""
mutable struct RealParam{T<:Real} <: AbstractParam{T}
    name::Symbol  # Name of the parameter

    val::T  # Current parameter value
    min_val::T  # Minimum parameter value
    max_val::T  # Maximum parameter value
    def_val::T  # Default value

    RealParam(name::Symbol, vdef::T, vmin::T, vmax::T) where{T<:Real} =
        new{T}(name, vdef, vmin, vmax, vdef)
    
end

function Base.copy(p::RealParam{T}) where{T<:Real}
    p_ = RealParam(p.name, p.def_val, p.min_val, p.max_val)
    p_.val = p.val
    return p_
end

const IntParam = RealParam{Int}
const FloatParam = RealParam{Float64}

get_param_name(par::RealParam) = par.name

get_param_type(par::RealParam{T}) where T = T

get_param_value(par::RealParam{T}) where T = par.val

set_param_default!(par::RealParam) = (par.val = par.def_val)

"""
    set_param_value!(p, v)

Set value of parameter `p` to `v`. Raises an error if `v` is not an admissible value.
"""
function set_param_value!(p::RealParam, v::T) where{T<:Real}

    if test_param_value(p, v)
        p.val = v
    else
        error("$(p.name) must be between $(p.min_val) and $(p.max_val).")
    end

    return nothing
end

"""
    test_param_value(p, v)

Return whether `v` is an admissible value for parameter `p` 
"""
test_param_value(p::RealParam, v::T) where{T<:Real} = (p.min_val <= v <= p.max_val)