@testset "CHOLMOD" begin

    A = SparseMatrixCSC{Float64, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    @testset "LDL" begin
        kkt = KKT.setup(A, KKT.K2(), KKT.TlpCholmod.Backend())
        KKT.run_ls_tests(A, kkt)
    end

    @testset "Cholesky" begin
        kkt = KKT.setup(A, KKT.K1(), KKT.TlpCholmod.Backend())
        KKT.run_ls_tests(A, kkt)
    end
end
