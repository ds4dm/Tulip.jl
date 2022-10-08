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
"""
    IPMData(pb::ProblemData, options::MatrixOptions)

Extract problem data to standard form.
"""
function IPMData(pb::ProblemData{T}, mfact::Factory) where{T}

    # Problem size
    m, n = pb.ncon, pb.nvar

    # Extract right-hand side and slack variables
    nzA = 0          # Number of non-zeros in A
    b = zeros(T, m)  # RHS
    sind = Int[]     # Slack row index
    sval = T[]       # Slack coefficient
    lslack = T[]     # Slack lower bound
    uslack = T[]     # Slack upper bound

    for (i, (lb, ub)) in enumerate(zip(pb.lcon, pb.ucon))
        if lb == ub
            # Equality row
            b[i] = lb

        elseif -T(Inf) == lb && T(Inf) == ub
            # Free row
            push!(sind, i)
            push!(sval, one(T))
            push!(lslack, -T(Inf))
            push!(uslack, T(Inf))
            b[i] = zero(T)

        elseif -T(Inf) == lb && isfinite(ub)
            # a'x <= b --> a'x + s = b
            push!(sind, i)
            push!(sval, one(T))
            push!(lslack, zero(T))
            push!(uslack, T(Inf))
            b[i] = ub

        elseif isfinite(lb) && ub == Inf
            # a'x >= b --> a'x - s = b
            push!(sind, i)
            push!(sval, -one(T))
            push!(lslack, zero(T))
            push!(uslack, T(Inf))
            b[i] = lb

        elseif isfinite(lb) && isfinite(ub)
            # lb <= a'x <= ub
            # Two options:
            # --> a'x + s = ub, 0 <= s <= ub - lb
            # --> a'x - s = lb, 0 <= s <= ub - lb 
            push!(sind, i)
            push!(sval, one(T))
            push!(lslack, zero(T))
            push!(uslack, ub - lb)
            b[i] = ub

        else
            error("Invalid bounds for row $i: [$lb, $ub]")
        end

        # This line assumes that there are no dupplicate coefficients in Arows
        # Numerical zeros will also be counted as non-zeros
        nzA += length(pb.arows[i].nzind)
    end

    nslack = length(sind)

    # Objective
    c = [pb.obj; zeros(T, nslack)]
    c0 = pb.obj0
    if !pb.objsense
        # Flip objective for maximization problem
        c .= .-c
        c0 = -c0
    end

    # Instantiate A
    aI = Vector{Int}(undef, nzA + nslack)
    aJ = Vector{Int}(undef, nzA + nslack)
    aV = Vector{T}(undef, nzA + nslack)

    # populate non-zero coefficients by column
    nz_ = 0
    for (j, col) in enumerate(pb.acols)
        for (i, aij) in zip(col.nzind, col.nzval)
            nz_ += 1

            aI[nz_] = i
            aJ[nz_] = j
            aV[nz_] = aij
        end
    end
    # populate slack coefficients
    for (j, (i, a)) in enumerate(zip(sind, sval))
        nz_ += 1
        aI[nz_] = i
        aJ[nz_] = n + j
        aV[nz_] = a
    end

    # At this point, we should have nz_ == nzA + nslack
    # If not, this means the data between rows and columns in `pb`
    # do not match each other
    nz_ == (nzA + nslack) || error("Found $(nz_) non-zero coeffs (expected $(nzA + nslack))")

    A = construct_matrix(mfact.T, m, n + nslack, aI, aJ, aV, mfact.options...)

    # Variable bounds
    l = [pb.lvar; lslack]
    u = [pb.uvar; uslack]
    
    return IPMData(A, b, pb.objsense, c, c0, l, u)
end
