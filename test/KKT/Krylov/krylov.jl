using Krylov

@testset "Krylov" begin
    for T in TvTYPES
        @testset "$T (CG)" begin
            A = SparseMatrixCSC{T, Int}([
                1 0 1 0;
                0 1 0 1
            ])
            kkt = KKT.setup(A, KKT.K1(), KKT.TlpKrylov.Backend(Krylov.CgSolver, Vector{T}))
            KKT.run_ls_tests(A, kkt)
        end
        @testset "$T (Minres)" begin
            A = SparseMatrixCSC{T, Int}([
                1 0 1 0;
                0 1 0 1
            ])
            kkt = KKT.setup(A, KKT.K1(), KKT.TlpKrylov.Backend(Krylov.MinresSolver, Vector{T}))
            KKT.run_ls_tests(A, kkt)
        end
    end
end
