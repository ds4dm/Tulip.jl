module Tulip

using LinearAlgebra
using Logging
using MutableArithmetics
using Printf
using SparseArrays
using TOML

using TimerOutputs

# Read Tulip version from Project.toml file
const TULIP_VERSION = let project = joinpath(@__DIR__, "..", "Project.toml")
    Base.include_dependency(project)
    VersionNumber(TOML.parsefile(project)["version"])
end
version() = TULIP_VERSION

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
