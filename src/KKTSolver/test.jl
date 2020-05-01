using Test
using LinearAlgebra

"""
    run_ls_tests(A, ls; atol)


"""
function run_ls_tests(
    A::AbstractMatrix{Tv},
    ls::AbstractKKTSolver{Tv};
    atol::Tv=sqrt(eps(Tv))
) where{Tv<:Real}

    # Check that required methods are implemented
    @testset "Required methods" begin
        Tls = typeof(ls)
        V = Vector{Tv}
        @test hasmethod(update_linear_solver!, Tuple{Tls, V, V, V})
        @test hasmethod(solve_augmented_system!, Tuple{V, V, Tls, V, V})
    end

    m, n = size(A)

    # Factorization/pre-conditionner update
    θ = ones(Tv, n)
    regP = ones(Tv, n)
    regD = ones(Tv, m)
    update_linear_solver!(ls, θ, regP, regD)

    # Solve linear system
    ξp = ones(Tv, m)
    ξd = ones(Tv, n)
    dx = zeros(Tv, n)
    dy = zeros(Tv, m)
    solve_augmented_system!(dx, dy, ls, ξp, ξd)

    # Check residuals
    rp = A * dx + regD .* dy - ξp
    rd = -dx .*(θ + regP) + A' * dy - ξd
    @testset "Residuals" begin
        @test norm(rp, Inf) <= atol
        @test norm(rd, Inf) <= atol
    end

    return nothing
end