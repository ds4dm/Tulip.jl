using Test

using LinearAlgebra
using SparseArrays

using Tulip
TLP = Tulip

const TvTYPES = [Float64, BigFloat]

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

@testset "Tulip Tests" begin 
    for f in test_files
        tp = joinpath(testdir, "$(f).jl")
        include(tp)
    end
end  # testset

# @testset "MathOptInterface Tests" begin
#     include("MOI_wrapper.jl")
# end