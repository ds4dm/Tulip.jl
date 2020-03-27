using Test

using LinearAlgebra
using SparseArrays

using Tulip
TLP = Tulip

const TvTYPES = [Float64, BigFloat]

# write your own tests here
const testdir = dirname(@__FILE__)

@testset "Tulip tests" begin
    @testset "Unit tests" begin
        include("Core/problemData.jl")
        include("Solvers/HSD.jl")
    end

    @testset "Linear Algebra" begin
        include("LinearAlgebra/linearAlgebra.jl")
    end

    @testset "Interfaces" begin
        # include("Interfaces/MOI_wrapper.jl")
    end

    @testset "Examples" begin
        include("examples.jl")
    end
end