import Krylov

@testset "Krylov" begin
    for T in TvTYPES
        @testset "$T" begin
            # Kyrlov solvers for symmetric positive-definite systems
            for f in [Krylov.cg, Krylov.minres, Krylov.minres_qlp]
                @testset "PD:  $f" begin
                    A = SparseMatrixCSC{T, Int}([
                        1 0 1 0;
                        0 1 0 1
                    ])
                    kkt = KKT.KrylovSPD(f, A)
                    KKT.run_ls_tests(A, kkt)
                end
            end
            # Kyrlov solvers for symmetric indefinite systems
            for f in [Krylov.minres, Krylov.minres_qlp]
                @testset "SID: $f" begin
                    A = SparseMatrixCSC{T, Int}([
                        1 0 1 0;
                        0 1 0 1
                    ])
                    kkt = KKT.KrylovSID(f, A)
                    KKT.run_ls_tests(A, kkt)
                end
            end
            # Kyrlov solvers for symmetric quasi-definite systems
            for f in [Krylov.tricg, Krylov.trimr]
                @testset "SQD: $f" begin
                    A = SparseMatrixCSC{T, Int}([
                        1 0 1 0;
                        0 1 0 1
                    ])
                    kkt = KKT.KrylovSQD(f, A)
                    KKT.run_ls_tests(A, kkt)
                end
            end
        end
    end
end
