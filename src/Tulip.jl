module Tulip

using LinearAlgebra
using Logging
using Printf
using SparseArrays

using TimerOutputs

# Linear algebra
include("LinearAlgebra/LinearAlgebra.jl")
import .TLPLinearAlgebra.construct_matrix
const TLA = TLPLinearAlgebra

# KKT solvers
include("KKT/KKT.jl")
using .KKT

# Commons data structures
# TODO: put this in a module
include("utils.jl")
include("status.jl")    # Termination and solution statuses
include("problemData.jl")
include("parameters.jl")
include("solution.jl")
include("attributes.jl")

# Presolve module
include("Presolve/Presolve.jl")

# IPM solvers
include("./IPM/IPM.jl")

# Model
include("model.jl")

# Interfaces
include("Interfaces/tulip_julia_api.jl")
include("Interfaces/MOI/MOI_wrapper.jl")

end # module