@testset "LAPACK" begin
    for Tv in TvTYPES
        @testset "$Tv" begin
            A = Matrix{Tv}([
                1 0 1 0;
                0 1 0 1
            ])
            ls = KKT.AbstractKKTSolver(KKT.Lapack(), KKT.NormalEquations(), A)
            KKT.run_ls_tests(A, ls)
        end
    end
end