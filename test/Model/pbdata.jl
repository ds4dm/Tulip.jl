# TODO: parametrize and run tests for mutiple types
function run_tests_pbdata(::Tv) where{Tv<:Real}

    @testset "Constructor" begin
        pb = TLP.ProblemData{Tv}()

        @test pb.constr_cnt == 0
        @test pb.var_cnt == 0
    end

    @testset "Variables" begin
        pb = TLP.ProblemData{Tv}()
        
        # Add variable
        vidx = TLP.new_variable_index!(pb)
        @test vidx.uuid == 1

        vdat = TLP.VarData{Tv}("x", zero(Tv), TLP.TLP_BND_LO, zero(Tv), Tv(Inf))
        v = TLP.Variable(vidx, vdat)
        TLP.add_variable!(pb, v)
        @test haskey(pb.vars, vidx)
        @test haskey(pb.var2con, vidx)

        # Add that variable again. Should raise an error
        @test_throws ErrorException TLP.add_variable!(pb, v)
    end

    @testset "Constraints" begin
        pb = TLP.ProblemData{Tv}()

        # Create new constraint
        cidx = TLP.new_constraint_index!(pb)
        @test cidx.uuid == 1

        cdat = TLP.LinConstrData{Tv}("c1", TLP.TLP_BND_RG, zero(Tv), oneunit(Tv))
        c = TLP.LinearConstraint(cidx, cdat)

        TLP.add_constraint!(pb, c)
        @test haskey(pb.constrs, cidx)
        @test haskey(pb.con2var, cidx)

        @test_throws ErrorException TLP.add_constraint!(pb, c)
    end

    @testset "Getters" begin
        
    end

end

@testset "ProblemData" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_pbdata(zero(Tv)) end
    end
end