const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Optimal" begin
    include(joinpath(examples_dir, "optimal.jl"))
    for T in TvTYPES
        @testset "$T" begin ex_optimal(T) end
    end
end
@testset "PrimalInfeas" begin
    include(joinpath(examples_dir, "infeasible.jl"))
    for T in TvTYPES
        @testset "$T" begin ex_infeasible(T) end
    end
end
@testset "DualInfeas" begin
    include(joinpath(examples_dir, "unbounded.jl"))
    for T in TvTYPES
        @testset "$T" begin ex_unbounded(T) end
    end
end
