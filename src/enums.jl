"""
    BoundType

"""
@enum BoundType begin
    TLP_BND_UP  # Upper-bounded
    TLP_BND_LO  # Lower-bounded
    TLP_BND_FX  # Fixed
    TLP_BND_FR  # Free
    TLP_BND_RG  # Range
end


function _check_bounds(bt::BoundType, lb::Tv, ub::Tv) where{Tv<:Real}
    if bt == TLP_BND_FX
        return Tv(-Inf)  < lb == ub  < Tv(Inf)
    elseif bt == TLP_BND_RG
        return Tv(-Inf)  < lb <= ub  < Tv(Inf)
    elseif bt == TLP_BND_UP
        return Tv(-Inf) == lb  < ub  < Tv(Inf)
    elseif bt == TLP_BND_LO
        return Tv(-Inf)  < lb  < ub == Tv(Inf)
    elseif bt == TLP_BND_FR
        return (Tv(-Inf) == lb) && (ub == Tv(Inf))
    else
        error("Unsupported bound type: $bt.")
    end

    return false
end