"""
    make_step!()

Compute next IP iterate for the HSD formulation.
"""
function make_step!(
    ::Val{1},
    model, env,
    A, F, b, c, uind, uval,
    θ, θwz,
    x, w, y, s, z, t, k
)
    # Affine-scaling direction

end

"""
    solve_newton_hsd!(...)

Solve Newton system
     A*dx           -b*dt                                   = rp
     U*dx   +  dw   -u*dt                                   = ru
                    -c*dt   +A'dy   +  ds   +U'dz           = rd
    -c'dx                   +b'dy           -u'dz   -  dk   = rg
     S*dx                           +X*ds                   = rxs
            +Z*dw                          +W*dz            = rwz
                     k*dt                           +t*dk   = rtk

This version assumes that the augmented system below has been solved first:

# Arguments
- A: Matrix
- `F`: Factorization of the normal equations matrix
- `b`: Right-hand side of primal constraints
- `c`: Primal objective term
- `uind`: Indices of upper-bounded variables
- `uval`: Upper bound values
- `θ`: Diagonal scaling term
- `θwz`: Diagonal scaling term
- `p, q, r, ρ`: Terms obtained in the preliminary augmented system solve
- `x, w, y, s, z, t, k`: Primal and dual iterates
- `dx, dw, dy, ds, dz, dt, dk`: Primal and dual directions, modified in-place
- `rp, rd, ru, rg, rxs, rwz, rtk`: 
"""
function solve_newton_hsd!(
    A, F, b, c, uind, uval, θ, θwz,
    p, q, r, ρ,
    x, w, y, s, z, t, k,
    dx, dw, dy, ds, dz, dt, dk,
    rp, ru, rd, rg::Float64, rxs, rwz, rtk::Float64
)
    # Solve reduced newton system
    solve_augsys_hsd!(
        A, F, θ, θwz, uind,
        dx, dy, dz,
        rp, 
        rd .- (rxs ./ x),
        ru .- (rwz ./ z)
    )

    # Compute Δτ
    dt.x = (
        rg + (rtk / t.x) + dot(c, dx) - dot(b, dy) + dot(uval, dz)
    ) / ρ
    
    dx .+= dt.x .* p
    dy .+= dt.x .* q
    dz .+= dt.x .* r

    # Compute Δs, Δκ, Δw
    # ds  .= (rxs - s  .* dx)   ./ x
    # dw  .= (rwz - w  .* dz)   ./ z
    ds  .= rxs
    ds .-= s .* dx
    ds ./= x

    dw  .= rwz
    dw .-= w .* dz
    dw ./= z

    dk.x = (rtk - k.x * dt.x)  / t.x

    return nothing
end

"""
    solve_augsys_hsd!(...)

Solve the augmented system below, and overwrite `dx, dy, dz` with the result.
    -(S/X)*dx  + A'dy  -    U'dz    = rd
         A*dx                       = rp
         U*dx          -(Z/W)*dz    = ru
The augmented system is solved by direct resolution of the normal equations.

# Arguments
- `A`: Matrix
- `F`: Factorization of the matrix `A*Θ*A'`
- `θ`
- `θwz`
- `uind`: Vector of indices of upper-bounded variables.
- `dx, dy, dz`: Vectors of unknowns, modified in-place.
- `rp, rd, ru`: Right-hand-side vectors.
"""
function solve_augsys_hsd!(
    A, F, θ, θwz, uind,
    dx, dy, dz,
    rp, rd, ru
)
    ru_ = ru .* θwz
    rp_ = copy(rd)
    @views rp_[uind] .-= ru_
    rp_ .*= θ
    mul!(dy, A, rp_)
    dy .+= rp

    dy .= (F \ dy)
    mul!(dx, transpose(A), dy)
    dx .-= rd
    @views dx[uind] .+= ru_
    dx .*= θ

    dz .= zero(eltype(dz))
    dz .-= ru
    @views dz .+= dx[uind]
    dz .*= θwz

    return nothing
end