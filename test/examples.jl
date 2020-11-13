const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Optimal" begin
    include(joinpath(examples_dir, "optimal.jl"))
    for T in TvTYPES
        @testset "$T" begin
            ex_optimal(T; BarrierAlgorithm=1, OutputLevel=0)
            ex_optimal(T; BarrierAlgorithm=2, OutputLevel=0)
        end
    end
end
@testset "Free vars" begin
    include(joinpath(examples_dir, "freevars.jl"))
    for T in TvTYPES
        @testset "$T" begin
            ex_freevars(T; BarrierAlgorithm=1, OutputLevel=0)
            ex_freevars(T; BarrierAlgorithm=2, OutputLevel=0)
        end
    end
end
@testset "PrimalInfeas" begin
    include(joinpath(examples_dir, "infeasible.jl"))
    for T in TvTYPES
        @testset "$T" begin ex_infeasible(T, OutputLevel=0) end
    end
end
@testset "DualInfeas" begin
    include(joinpath(examples_dir, "unbounded.jl"))
    for T in TvTYPES
        @testset "$T" begin ex_unbounded(T, OutputLevel=0) end
    end
end
