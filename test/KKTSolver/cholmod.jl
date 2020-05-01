@testset "CHOLMOD" begin

    A = SparseMatrixCSC{Float64, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    @testset "LDL" begin
        ls = KKT.AbstractLinearSolver(KKT.Cholmod(), KKT.AugmentedSystem(), A)
        KKT.run_ls_tests(A, ls)
    end
    
    @testset "Cholesky" begin
        ls = KKT.AbstractLinearSolver(KKT.Cholmod(), KKT.NormalEquations(), A)
        KKT.run_ls_tests(A, ls)
    end

end