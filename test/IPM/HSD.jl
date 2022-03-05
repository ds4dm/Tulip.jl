function run_tests_hsd(T::Type)

    Tv = Vector{T}

    params = TLP.IPMOptions{T}()
    kkt_options = TLP.KKTOptions{T}()

    @testset "step length" begin
        m, n, p = 2, 2, 1
        pt = TLP.Point{T, Tv}(m, n, p, hflag=true)
        pt.x  .= one(T)
        pt.xl .= one(T)
        pt.xu .= one(T)
        pt.y  .= zero(T)
        pt.zl .= zero(T)
        pt.zu .= zero(T)
        pt.τ   = one(T)
        pt.κ   = one(T)
        pt.μ   = one(T)

        d = TLP.Point{T, Tv}(m, n, p, hflag=true)
        d.x  .= one(T)
        d.xl .= one(T)
        d.xu .= one(T)
        d.y  .= zero(T)
        d.zl .= zero(T)
        d.zu .= zero(T)
        d.τ   = one(T)
        d.κ   = one(T)
        d.μ   = one(T)

        # Max step length for a single (x, d)
        @inferred TLP.max_step_length(pt.x, d.x)
        @test TLP.max_step_length(ones(T, 1), ones(T, 1)) == T(Inf)
        @test TLP.max_step_length(ones(T, 1), -ones(T, 1)) ≈ one(T)
        @test TLP.max_step_length(zeros(T, 1), -ones(T, 1)) ≈ zero(T)
        @test TLP.max_step_length(zeros(T, 1), ones(T, 1)) == T(Inf)

        # Max step length for the whole primal-dual point
        @inferred TLP.max_step_length(pt, d)
        @test TLP.max_step_length(pt, d) ≈ one(T)
    end

    # Simple example:
    #=
        min     x1 - x2
        s.τ.    x1 + x2 = 1
                x1 - x2 = 0
                0 <= x1 <= 2
                0 <= x2 <= 2
    =#
    m, n = 2, 2
    p = 2 * n
    A = Matrix{T}([
        [1    1];
        [1   -1]
    ])
    b = Vector{T}([1,  0])
    c = Vector{T}([1, -1])
    c0 = zero(T)
    l = Vector{T}([0, 0])
    u = Vector{T}([2, 2])
    dat = Tulip.IPMData(A, b, true, c, c0, l, u)

    hsd = TLP.HSD(dat, kkt_options)

    # Primal-dual optimal solution
    # x1 = x2 = 0.5; xl = 0.5; xu = 1.5; τ = 1
    # y1 = 0, y2 = 1; zl = zu = 0; κ = 0
    hsd.pt.x  .= T.([1 // 2, 1 // 2])
    hsd.pt.xl .= T.([1 // 2, 1 // 2])
    hsd.pt.xu .= T.([3 // 2, 3 // 2])
    hsd.pt.y  .= T.([0, 1])
    hsd.pt.zl .= T.([0, 0])
    hsd.pt.zu .= T.([0, 0])

    hsd.pt.τ  = 1
    hsd.pt.κ  = 0
    hsd.pt.μ  = 0

    ϵ = sqrt(eps(T))
    TLP.compute_residuals!(hsd)

    @testset "Convergence" begin
        hsd.solver_status = TLP.Trm_Unknown
        TLP.update_solver_status!(hsd, ϵ, ϵ, ϵ, ϵ)
        @test hsd.solver_status == TLP.Trm_Optimal

        # TODO: dual infeasible

        # TODO: primal infeasible

        # TODO: ill-posed

    end
end

function test_hsd_residuals(T::Type)
    # Simple example:
    #=
        min     x1 - x2
        s.τ.    x1 + x2 = 1
                x1 - x2 = 0
                0 <= x1 <= 2
                0 <= x2 <= 2
    =#
    kkt_options = TLP.KKTOptions{T}()
    A = Matrix{T}([
        [1    1];
        [1   -1]
    ])
    b = Vector{T}([1,  0])
    c = Vector{T}([1, -1])
    c0 = zero(T)
    l = Vector{T}([0, 0])
    u = Vector{T}([2, 2])
    dat = Tulip.IPMData(A, b, true, c, c0, l, u)

    hsd = TLP.HSD(dat, kkt_options)
    pt = hsd.pt
    res = hsd.res

    # Primal-dual solution
    x  = pt.x  .= T[3, 5]
    xl = pt.xl .= T[1, 8]
    xu = pt.xu .= T[2, 1]
    y  = pt.y  .= T[10, -2]
    zl = pt.zl .= T[2, 1]
    zu = pt.zu .= T[5, 7]
    τ  = pt.τ  = T(1//2)
    κ  = pt.κ  = T(1//10)
    μ  = pt.μ  = 0

    TLP.compute_residuals!(hsd)

    @test res.rp ≈ (τ .* b) - A * x
    @test res.rl ≈ (τ .* l) - (x - xl)
    @test res.ru ≈ (τ .* u) - (x + xu)
    @test res.rd ≈ (τ .* c) - A' * y - zl + zu
    @test res.rg ≈ c'x - (b'y + l'zl - u'zu) + κ

    @test res.rp_nrm == norm(res.rp, Inf)
    @test res.rl_nrm == norm(res.rl, Inf)
    @test res.ru_nrm == norm(res.ru, Inf)
    @test res.rd_nrm == norm(res.rd, Inf)
    @test res.rg_nrm == norm(res.rg, Inf)

    return nothing
end

@testset "HSD" begin
    @testset "$T" for T in TvTYPES
        run_tests_hsd(T)
        test_hsd_residuals(T)
    end
end
