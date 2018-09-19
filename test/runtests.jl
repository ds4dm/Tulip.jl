using Tulip
using Test

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "params",
    "env",
    "LinearAlgebra/denseBlockAngular",
    "LinearAlgebra/linearAlgebra",
    "model",
    "mathprogbase",
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end