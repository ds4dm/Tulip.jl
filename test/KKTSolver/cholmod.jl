@testset "CHOLMOD" begin

    A = SparseMatrixCSC{Float64, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    @testset "LDL" begin
        ls = KKT.CholmodSolver(A, normal_equations=false)
        KKT.run_ls_tests(A, ls)
    end
    
    @testset "Cholesky" begin
        ls = KKT.CholmodSolver(A, normal_equations=true)
        KKT.run_ls_tests(A, ls)
    end

end