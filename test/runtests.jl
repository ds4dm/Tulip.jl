using LinearAlgebra
using SparseArrays
using Test
using TOML

using Tulip
TLP = Tulip

const TvTYPES = [Float32, Float64, BigFloat]

# Check That Tulip.version() matches what's in the Project.toml
tlp_ver = Tulip.version()
toml_ver = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
@test tlp_ver == toml_ver

@testset "Tulip" begin

@testset "Unit tests" begin

    @testset "Core" begin
        include("Core/problemData.jl")
    end

    @testset "IPM" begin
        include("IPM/HSD.jl")
        include("IPM/MPC.jl")
    end

    @testset "KKT" begin
        include("KKT/KKT.jl")
    end

    @testset "Presolve" begin
        include("Presolve/presolve.jl")
    end
end  # UnitTest

@testset "Examples" begin
    include("examples.jl")
end

@testset "Interfaces" begin
    include("Interfaces/julia_api.jl")
end

@testset "MOI" begin
    include("Interfaces/MOI_wrapper.jl")
end

# @testset "Convex Problem Depot tests" begin
#     for T in TvTYPES
#         @testset "$T" begin
#             Convex.ProblemDepot.run_tests(;  exclude=[r"mip", r"exp", r"socp", r"sdp"], T = T) do problem
#                 Convex.solve!(problem, () -> Tulip.Optimizer{T}())
#             end
#         end
#     end
# end

end  # Tulip tests
