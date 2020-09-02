import Krylov

@testset "Krylov" begin
    for T in TvTYPES
        @testset "$T" begin
            for f in [Krylov.cg, Krylov.minres, Krylov.minres_qlp]
                @testset "$f" begin
                    A = SparseMatrixCSC{T, Int}([
                        1 0 1 0;
                        0 1 0 1
                    ])
                    kkt = KKT.KrylovPDSolver(f, A)
                    KKT.run_ls_tests(A, kkt)
                end
            end
            for f in [Krylov.tricg, Krylov.trimr]
                @testset "$f" begin
                    A = SparseMatrixCSC{T, Int}([
                        1 0 1 0;
                        0 1 0 1
                    ])
                    kkt = KKT.KrylovSQDSolver(f, A)
                    KKT.run_ls_tests(A, kkt)
                end
            end
        end
    end
end