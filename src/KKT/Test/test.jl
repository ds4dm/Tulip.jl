using Test
using LinearAlgebra

"""
    run_ls_tests(A, kkt; atol)


"""
function run_ls_tests(
    A::AbstractMatrix{T},
    kkt::AbstractKKTSolver{T};
    atol::T=sqrt(eps(T))
) where{T}

    # Check that required methods are implemented
    @testset "Required methods" begin
        Tls = typeof(kkt)
        V = Vector{T}
        @test hasmethod(update!, Tuple{Tls, V, V, V})
        @test hasmethod(solve!, Tuple{V, V, Tls, V, V})
    end

    m, n = size(A)

    # Factorization/pre-conditionner update
    θ = ones(T, n)
    regP = ones(T, n)
    regD = ones(T, m)
    update!(kkt, θ, regP, regD)

    # Solve linear system
    ξp = ones(T, m)
    ξd = ones(T, n)
    dx = zeros(T, n)
    dy = zeros(T, m)
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