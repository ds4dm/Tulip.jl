function test_krylov_sqd(T, ksolver)
    A = SparseMatrixCSC{T, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    kkt = KKT.setup(A, KKT.K2(), KKT.TlpKrylov.Backend(ksolver, Vector{T}))
    KKT.run_ls_tests(A, kkt)

    return nothing
end

@testset "SQD" begin
    for T in TvTYPES, ksolver in [TricgWorkspace, TrimrWorkspace]
        @testset "$ksolver ($T)" begin
            test_krylov_sqd(T, ksolver)
        end
    end
end
