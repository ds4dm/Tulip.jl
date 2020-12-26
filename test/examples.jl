const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Optimal" begin
    include(joinpath(examples_dir, "optimal.jl"))
    for T in TvTYPES
        @testset "$T" begin
            ex_optimal(T; OutputLevel=1, IPM_Factory=Tulip.Factory(Tulip.HSD))
            ex_optimal(T; OutputLevel=1, IPM_Factory=Tulip.Factory(Tulip.MPC))
        end
    end
end
@testset "Free vars" begin
    include(joinpath(examples_dir, "freevars.jl"))
    for T in TvTYPES
        @testset "$T" begin
            ex_freevars(T; OutputLevel=1, IPM_Factory=Tulip.Factory(Tulip.HSD))
            ex_freevars(T; OutputLevel=1, IPM_Factory=Tulip.Factory(Tulip.MPC))
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

@testset "Optimal Float32" begin
    include(joinpath(examples_dir, "optimal_other_type.jl"))
end
