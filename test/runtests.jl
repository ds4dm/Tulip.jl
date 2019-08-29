using Test

using LinearAlgebra
using SparseArrays

using Tulip
TLP = Tulip

const TvTYPES = [Float32, Float64, BigFloat, Rational]

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    "env",
    "bounds",
    "LinearAlgebra/linearAlgebra",
    "Solvers/HSD",
    "Model/model",
    "reader",
    "examples"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end