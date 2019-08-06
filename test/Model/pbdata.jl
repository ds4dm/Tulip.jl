function run_tests_pbdata()

    @testset "Constructor" begin
        pb = TLP.ProblemData{Float64}()

        @test pb.constr_cnt == 0
        @test pb.var_cnt == 0
    end

    @testset "Variables" begin
        pb = TLP.ProblemData{Float64}()
        
        # Add variable
        vidx = TLP.new_variable_index!(pb)
        @test vidx.uuid == 1
        v = TLP.Variable{Float64}(vidx)
        TLP.add_variable!(pb, v)
        @test haskey(pb.vars, vidx)
        @test haskey(pb.var2con, vidx)

        # Add that variable again. Should raise an error
        @test_throws ErrorException TLP.add_variable!(pb, v)
    end

    @testset "Constraints" begin
        pb = TLP.ProblemData{Float64}()

        # Create new constraint
        cidx = TLP.new_constraint_index!(pb)
        @test cidx.uuid == 1

        c = TLP.LinearConstraint{Float64}(cidx)
        TLP.add_constraint!(pb, c)
        @test haskey(pb.constrs, cidx)
        @test haskey(pb.con2var, cidx)

        @test_throws ErrorException TLP.add_constraint!(pb, c)
    end

    @testset "Getters" begin
        
    end

end

@testset "ProblemData" begin
    run_tests_pbdata()
end