import Tulip: HSDSolver, optimize!

using Random

function run_test_hsd(::Tv) where{Tv<:Real}
    
end

@testset "HSD" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_test_hsd(zero(Tv)) end
    end
end