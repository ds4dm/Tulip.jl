"""
    BoundType

* `TLP_UP`: Upper-bounded variable ``-\\infty \\leq x \\leq u``.
* `TLP_LO`: Lower-bounded variable ``l \\leq x \\leq \\infty``.
* `TLP_FX`: Fixed variable ``l = x = u``.
* `TLP_FR`: Free variable ``-\\infty \\leq x \\leq \\infty``.
* `TLP_RG`: Range variable ``l \\leq x \\leq u``.
"""
@enum BoundType begin
    TLP_UP  # Upper-bounded
    TLP_LO  # Lower-bounded
    TLP_FX  # Fixed
    TLP_FR  # Free
    TLP_RG  # Range
end


function bound_type(lb::Tv, ub::Tv)::BoundType where{Tv<:Real}

    if Tv(-Inf)  < lb == ub  < Tv(Inf)
        return TLP_FX

    elseif Tv(-Inf)  < lb <= ub  < Tv(Inf)
        return TLP_RG

    elseif Tv(-Inf) == lb  < ub  < Tv(Inf)
        return TLP_UP

    elseif Tv(-Inf)  < lb  < ub == Tv(Inf)
        return TLP_LO

    elseif (Tv(-Inf) == lb) && (ub == Tv(Inf))
        return TLP_FR

    else
        error("Invalid bounds: [$lb, $ub].")
    end
end