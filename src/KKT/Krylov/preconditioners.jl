abstract type AbstractPreconditioner{T} end

update_preconditioner(::AbstractKKTSolver) = nothing

struct IdentityPreconditioner end

op(::IdentityPreconditioner) = LO.opEye()

include("jacobi.jl")