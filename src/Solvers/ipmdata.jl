"""
    IPMData{T, Tv, Ta}

Holds data about an interior point method.

The problem is represented as
```
min   c'x + c0
s.t.  A x = b
      l ≤ x ≤ u
```
where `l`, `u` may take infinite values.
"""
struct IPMData{T, Tv, Tb, Ta}

    # Problem size
    nrow::Int
    ncol::Int

    # Objective
    objsense::Bool  # min (true) or max (false)
    c0::T
    c::Tv

    # Constraint matrix
    A::Ta

    # RHS
    b::Tv

    # Variable bounds (may contain infinite values)
    l::Tv
    u::Tv
    # Variable bound flags (we template with `Tb` to ease GPU support)
    # These should be vectors of the same type as `l`, `u`, but `Bool` eltype.
    # They should not be passed as arguments, but computed at instantiation as
    # `lflag = isfinite.(l)` and `uflag = isfinite.(u)`
    lflag::Tb
    uflag::Tb

    function IPMData(
        A::Ta, b::Tv, objsense::Bool, c::Tv, c0::T, l::Tv, u::Tv
    ) where{T, Tv<:AbstractVector{T}, Ta<:AbstractMatrix{T}}
        nrow, ncol = size(A)

        lflag = isfinite.(l)
        uflag = isfinite.(u)
        Tb = typeof(lflag)

        return new{T, Tv, Tb, Ta}(
            nrow, ncol,
            objsense, c0, c,
            A, b, l, u, lflag, uflag
        )
    end
end

# TODO: extract IPM data from presolved problem