@testset "CHOLMOD" begin

    A = SparseMatrixCSC{Float64, Int}([
        1 0 1 0;
        0 1 0 1
    ])

    @testset "LDL" begin
        kkt = KKT.CholmodSolver(A, normal_equations=false)
        KKT.run_ls_tests(A, kkt)
    end
    
    @testset "Cholesky" begin
        kkt = KKT.CholmodSolver(A, normal_equations=true, variant=false)
        KKT.run_ls_tests(A, kkt)
    end

    @testset "Cholesky" begin
        kkt = KKT.CholmodSolver(A, normal_equations=true, variant=true)
        KKT.run_ls_tests(A, kkt)
    end

end