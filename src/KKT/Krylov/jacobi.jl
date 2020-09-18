struct JacobiPreconditioner{T, Tv} <: AbstractPreconditioner{T}
    d::Tv

    JacobiPreconditioner(d::Tv) where{T, Tv<:AbstractVector{T}} = new{T, Tv}(d)
end

op(P::JacobiPreconditioner) = Diagonal(P.d)

# Code for Jacobi preconditioners
function update_preconditioner(kkt::KrylovSPD{T, F, Tv, Ta, Tp}) where{T, F, Tv, Ta<:SparseMatrixCSC, Tp<:JacobiPreconditioner{T, Tv}}
    
    A = kkt.A
    d = kkt.P.d

    # S = A (Θ⁻¹ + Rp)⁻¹ A' + Rd, so its diagonal writes
    # S[i, i] = Σ_j=1..n D[j, j] * A[i, j] ^2,
    # where D = (Θ⁻¹ + Rp)⁻¹

    _d = one(T) ./ (kkt.θ .+ kkt.regP)

    rows = rowvals(A)
    vals = nonzeros(A)
    m, n = size(A)

    d .= zero(T)

    for j = 1:n
        _a = zero(T)
        dj = _d[j]
        for i in nzrange(A, j)
            row = rows[i]
            a_ij = vals[i]

            d[row] += dj * a_ij ^2
        end
    end
    
    # I think putting this afterwards is more numerics-friendly
    d .+= kkt.regD
    
    # Invert
    d .\= one(T)

    return nothing
end