using Test

using LinearAlgebra
using SparseArrays

using Tulip
TLP = Tulip

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "Model/model",
    "params",
    "env",
    "LinearAlgebra/linearAlgebra",
    # "hsd",
    # "model",
    # "mathprogbase",
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end