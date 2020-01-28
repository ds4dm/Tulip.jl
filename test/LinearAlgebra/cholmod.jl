@testset "CHOLMOD" begin

    A = SparseMatrixCSC{Float64, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    @testset "LDL" begin
        ls = TLA.AbstractLinearSolver(TLA.Cholmod(), TLA.AugmentedSystem(), A)
        TLA.run_ls_tests(A, ls)
    end
    
    @testset "Cholesky" begin
        ls = TLA.AbstractLinearSolver(TLA.Cholmod(), TLA.NormalEquations(), A)
        TLA.run_ls_tests(A, ls)
    end

end