function test_krylov_sid(T, ksolver)
    A = SparseMatrixCSC{T, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    kkt = KKT.setup(A, KKT.K2(), KKT.TlpKrylov.Backend(ksolver, Vector{T}))
    KKT.run_ls_tests(A, kkt)

    return nothing
end

@testset "SID" begin
    for T in TvTYPES, ksolver in [MinresWorkspace, MinresQlpWorkspace, SymmlqWorkspace]
        @testset "$ksolver ($T)" begin
            test_krylov_sid(T, ksolver)
        end
    end
end
