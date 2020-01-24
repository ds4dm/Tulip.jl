import Tulip: HSDSolver, optimize!

function run_tests_hsd(::Tv) where{Tv<:Real}

    env = TLP.Env{Tv}()

    @testset "step length" begin
        m, n, p = 2, 2, 1
        pt = TLP.Point{Tv}(m, n, p)
        pt.x .= one(Tv)
        pt.w .= one(Tv)
        pt.t  = one(Tv)
        pt.s .= one(Tv)
        pt.z .= one(Tv)
        pt.k  = one(Tv)
        pt.μ  = one(Tv)

        d = TLP.Point{Tv}(m, n, p)
        d.x .= -one(Tv)
        d.w .= -one(Tv)
        d.t  = -one(Tv)
        d.s .= -one(Tv)
        d.z .= -one(Tv)
        d.k  = -one(Tv)
        d.μ  =  one(Tv)

        # Max step length for a single (x, d)
        @inferred TLP.max_step_length(pt.x, d.x)
        @test TLP.max_step_length(ones(Tv, 1), ones(Tv, 1)) == Tv(Inf)
        @test TLP.max_step_length(ones(Tv, 1), -ones(Tv, 1)) ≈ one(Tv)
        @test TLP.max_step_length(zeros(Tv, 1), -ones(Tv, 1)) ≈ zero(Tv)
        @test TLP.max_step_length(zeros(Tv, 1), ones(Tv, 1)) == Tv(Inf)

        # Max step length for the whole primal-dual point
        @inferred TLP.max_step_length(pt, d)
        @test TLP.max_step_length(pt, d) ≈ one(Tv)
    end

    @testset "Residuals" begin
        # Simple example:
        #=
            min     x1 - x2
            s.t.    x1 + x2 = 1
                    x1 - x2 = 0
                    0 <= x1 <= 2
                    0 <= x2 <= 2
        =#
        m, n, p = 2, 2, 2
        A = Matrix{Tv}([
            [1.0    1.0];
            [1.0   -1.0]
        ])
        b = Vector{Tv}([1.0, 0.0])
        c = Vector{Tv}([1.0, -1.0])
        c0 = zero(Tv)
        uind = [1, 2]
        uval = Vector{Tv}([2.0, 2.0])

        hsd = TLP.HSDSolver{Tv}(env, m, n, p, A, b, c, c0, uind, uval)

        # Primal-dual optimal solution
        # x1 = x2 = 0.5; w1 = w2 = 1.5; t = 1
        # y1 = 0, y2 = 1; s1 = s2 = 0; z1 = z2 = 0; k = 0
        pt = TLP.Point{Tv}(m, n, p)
        pt.x .= Tv.([0.5, 0.5])
        pt.w .= Tv.([1.5, 1.5])
        pt.t  = 1
        pt.y .= Tv.([0.0, 1.0])
        pt.s .= 0
        pt.z .= 0
        pt.k  = 0
        pt.μ  = 0

        res = TLP.Residuals{Tv}(
            zeros(Tv, m), zeros(Tv, p), zeros(Tv, n), zero(Tv),
            zero(Tv), zero(Tv), zero(Tv), zero(Tv)
        )

        @inferred TLP.compute_residuals!(hsd, res, pt, A, b, c, c0, uind, uval)
        TLP.compute_residuals!(hsd, res, pt, A, b, c, c0, uind, uval)
        
        @test res.rp_nrm ≈ zero(Tv)
        @test res.ru_nrm ≈ zero(Tv)
        @test res.rd_nrm ≈ zero(Tv)
        # @test res.rg_nrm == zero(Tv)
        
    end

    @testset "Convergence" begin

        # Optimal case
        #=
            min     x1 - x2
            s.t.    x1 + x2 = 1
                    x1 - x2 = 0
                    0 <= x1 <= 2
                    0 <= x2 <= 2
        =#
        m, n, p = 2, 2, 2
        A = Matrix{Tv}([
            [1.0    1.0];
            [1.0   -1.0]
        ])
        b = Vector{Tv}([1.0, 0.0])
        c = Vector{Tv}([1.0, -1.0])
        c0 = zero(Tv)
        uind = [1, 2]
        uval = Vector{Tv}([2.0, 2.0])

        hsd = TLP.HSDSolver{Tv}(env, m, n, p, A, b, c, c0, uind, uval)

        # Primal-dual optimal solution
        # x1 = x2 = 0.5; w1 = w2 = 1.5; t = 1
        # y1 = 0, y2 = 1; s1 = s2 = 0; z1 = z2 = 0; k = 0
        pt = TLP.Point{Tv}(m, n, p)
        pt.x .= Tv.([0.5, 0.5])
        pt.w .= Tv.([1.5, 1.5])
        pt.t  = 1
        pt.y .= Tv.([0.0, 1.0])
        pt.s .= 0
        pt.z .= 0
        pt.k  = 0
        pt.μ  = 0
        res = TLP.Residuals{Tv}(
            zeros(Tv, m), zeros(Tv, p), zeros(Tv, n), zero(Tv),
            zero(Tv), zero(Tv), zero(Tv), zero(Tv)
        )
        TLP.compute_residuals!(hsd, res, pt, A, b, c, c0, uind, uval)

        hsd.solver_status = TLP.Trm_Unknown
        TLP.update_solver_status!(hsd, pt, res, A, b, c, c0, uind, uval,
            Tv(1e-8), Tv(1e-8), Tv(1e-8), Tv(1e-8)
        )
        @test hsd.solver_status == TLP.Trm_Optimal

        # TODO: dual infeasible

        # TODO: primal infeasible

        # TODO: ill-posed

    end
end

@testset "HSD" begin
    for Tv in TvTYPES
        @testset "$Tv" begin run_tests_hsd(zero(Tv)) end
    end
end