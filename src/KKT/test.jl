using Test
using LinearAlgebra

"""
    run_ls_tests(A, kkt; atol)


"""
function run_ls_tests(
    A::AbstractMatrix{Tv},
    kkt::AbstractKKTSolver{Tv};
    atol::Tv=sqrt(eps(Tv))
) where{Tv<:Real}

    # Check that required methods are implemented
    @testset "Required methods" begin
        Tls = typeof(kkt)
        V = Vector{Tv}
        @test hasmethod(update!, Tuple{Tls, V, V, V})
        @test hasmethod(solve!, Tuple{V, V, Tls, V, V})
    end

    m, n = size(A)

    # Factorization/pre-conditionner update
    θ = ones(Tv, n)
    regP = ones(Tv, n)
    regD = ones(Tv, m)
    update!(kkt, θ, regP, regD)

    # Solve linear system
    ξp = ones(Tv, m)
    ξd = ones(Tv, n)
    dx = zeros(Tv, n)
    dy = zeros(Tv, m)
    solve!(dx, dy, kkt, ξp, ξd)

    # Check residuals
    rp = A * dx + regD .* dy - ξp
    rd = -dx .*(θ + regP) + A' * dy - ξd
    @testset "Residuals" begin
        @test norm(rp, Inf) <= atol
        @test norm(rd, Inf) <= atol
    end

    return nothing
end