const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Optimal" begin
    include(joinpath(examples_dir, "optimal.jl"))
    for Tv in TvTYPES
        @testset "$Tv" begin ex_optimal(Tv) end
    end
end
@testset "PrimalInfeas" begin
    include(joinpath(examples_dir, "infeasible.jl"))
    for Tv in TvTYPES
        @testset "$Tv" begin ex_infeasible(Tv) end
    end
end
@testset "DualInfeas" begin
    include(joinpath(examples_dir, "unbounded.jl"))
    for Tv in TvTYPES
        @testset "$Tv" begin ex_unbounded(Tv) end
    end
end
