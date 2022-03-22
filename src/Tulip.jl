module Tulip

using LinearAlgebra
using Logging
using Printf
using SparseArrays
using TOML

using TimerOutputs

const _TULIP_VERSION = Ref{VersionNumber}()

function __init__()
    # Read Tulip version from Project.toml file
    tlp_ver = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
    _TULIP_VERSION[] = tlp_ver
end

version() = _TULIP_VERSION[]

include("utils.jl")

# Linear algebra
include("LinearAlgebra/LinearAlgebra.jl")
import .TLPLinearAlgebra.construct_matrix
const TLA = TLPLinearAlgebra

# KKT solvers
include("KKT/KKT.jl")
using .KKT

# Commons data structures
# TODO: put this in a module

include("status.jl")    # Termination and solution statuses
include("problemData.jl")
include("solution.jl")
include("attributes.jl")

# Presolve module
include("Presolve/Presolve.jl")

# IPM solvers
include("./IPM/IPM.jl")


include("parameters.jl")

# Model
include("model.jl")

# Interfaces
include("Interfaces/tulip_julia_api.jl")
include("Interfaces/MOI/MOI_wrapper.jl")

end # module
