function test_krylov_spd(T, ksolver)
    A = SparseMatrixCSC{T, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    kkt = KKT.setup(A, KKT.K1(), KKT.TlpKrylov.Backend(ksolver, Vector{T}))
    KKT.run_ls_tests(A, kkt)

    return nothing
end

@testset "SPD" begin
    for T in TvTYPES, ksolver in [CgWorkspace, CarWorkspace]
        @testset "$ksolver ($T)" begin
            test_krylov_spd(T, ksolver)
        end
    end
end
