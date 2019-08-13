const examples_dir = joinpath(@__FILE__, "../../examples")

@testset "Examples" begin
    include(joinpath(examples_dir, "infeasible.jl"))
    include(joinpath(examples_dir, "unbounded.jl"))
end