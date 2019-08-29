const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Examples" begin
    @testset "Optimal"      begin include(joinpath(examples_dir, "optimal.jl")) end
    @testset "PrimalInfeas" begin include(joinpath(examples_dir, "infeasible.jl")) end
    @testset "DualInfeas"   begin include(joinpath(examples_dir, "unbounded.jl")) end
end