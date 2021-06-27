"""
    DefaultKKTSystem

Default KKT system setting. Currently equivalent to [K2](@ref)
"""
struct DefaultKKTSystem <: AbstractKKTSystem end

"""
    K2

Augmented system
```
    [-(Θ⁻¹ + Rp)   Aᵀ] [dx] = [ξd]
    [   A          Rd] [dy]   [ξp]
```
where `ξd` and `ξp` are given right-hand side.
"""
struct K2 <: AbstractKKTSystem end

"""
    K1

Normal equations system
```
    (A * ((Θ⁻¹ + Rp)⁻¹ * Aᵀ + Rd) dy = ξp + A * (θ⁻¹ + Rp)⁻¹ * ξd
                                  dx = (Θ⁻¹ + Rp)⁻¹ * (Aᵀ * dy - ξd)
```
"""
struct K1 <: AbstractKKTSystem end
