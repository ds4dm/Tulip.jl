"""
    BoundType

"""
@enum BoundType begin
    TLP_UP  # Upper-bounded
    TLP_LO  # Lower-bounded
    TLP_FX  # Fixed
    TLP_FR  # Free
    TLP_RG  # Range
end


function _check_bounds(bt::BoundType, lb::Tv, ub::Tv) where{Tv<:Real}
    if bt == TLP_FX
        return Tv(-Inf)  < lb == ub  < Tv(Inf)
    elseif bt == TLP_RG
        return Tv(-Inf)  < lb <= ub  < Tv(Inf)
    elseif bt == TLP_UP
        return Tv(-Inf) == lb  < ub  < Tv(Inf)
    elseif bt == TLP_LO
        return Tv(-Inf)  < lb  < ub == Tv(Inf)
    elseif bt == TLP_FR
        return (Tv(-Inf) == lb) && (ub == Tv(Inf))
    else
        error("Unsupported bound type: $bt.")
    end

    return false
end