@testset "LAPACK" begin
    for Tv in TvTYPES
        @testset "$Tv" begin
            A = Matrix{Tv}([
                1 0 1 0;
                0 1 0 1
            ])
            kkt = KKT.Dense_SymPosDef(A)
            KKT.run_ls_tests(A, kkt)
        end
    end
end