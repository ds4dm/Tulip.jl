function test_linalg_dense(Tv::Type{<:Real})
    
    A = Matrix{Tv}([
        1 0 1 0;
        0 1 0 1
    ])
    test_linalg(A)

    # TODO: test linear solvers more extensively (?)
    return nothing
end

@testset "Matrix" begin
    for Tv in TvTYPES
        @testset "$Tv" begin
            test_linalg_dense(Tv)
        end    
    end
end  # testset