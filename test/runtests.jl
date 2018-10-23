using Tulip
using Test

using LinearAlgebra
using Random
using SparseArrays

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "params",
    "env",
    "LinearAlgebra/denseBlockAngular",
    "LinearAlgebra/unitBlockAngular",
    "LinearAlgebra/linearAlgebra",
    "hsd",
    "model",
    "mathprogbase",
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end