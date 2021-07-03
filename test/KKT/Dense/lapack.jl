@testset "LAPACK" begin
    for T in TvTYPES
        @testset "$T" begin
            A = Matrix{T}([
                1 0 1 0;
                0 1 0 1
            ])
            kkt = KKT.setup(A, KKT.K1(), KKT.TlpDense.Backend())
            KKT.run_ls_tests(A, kkt)
        end
    end
end
