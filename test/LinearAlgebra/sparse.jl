function test_linalg_sparse(Tv::Type{<:Real})
    
    A = SparseMatrixCSC{Tv, Int}([
        1 0 1 0;
        0 1 0 1
    ])
    test_linalg(A)

    # TODO: test linear solvers more extensively (?)
    return nothing
end

@testset "SparseMatrixCSC" begin
    for Tv in TvTYPES
        Tv == Float64 || continue  # Only do Float64 for now
        @testset "$Tv" begin
            test_linalg_sparse(Tv)
        end    
    end
end  # testset