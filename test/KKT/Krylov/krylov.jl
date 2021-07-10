using Krylov

@testset "Krylov" begin
    include("spd.jl")
    include("sid.jl")
    include("sqd.jl")
end
