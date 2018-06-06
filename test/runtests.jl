using Tulip
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("../src/Tulip.jl")

# write your own tests here
const testdir = dirname(@__FILE__)

const test_files = [
    # include test file name here (without .jl extension)
    # "model",
    # "tst_basic",
    # "denseBlockAngular",
    # "linear_algebra",
    "MPBInterface"
]

for f in test_files
    tp = joinpath(testdir, "$(f).jl")
    include(tp)
end
