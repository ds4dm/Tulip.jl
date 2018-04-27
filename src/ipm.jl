"""
solve(A, b, c, x0, y0, s0)

Infeasible predictor-corrector Interior-Point algorithm.

Solves the (continuous) linear program
    ``min c^{T} x`` s.t. ``Ax = b``
using an infeasible Primal-Dual Interior-Point algorithm.

# Arguments
- `A`: constraint matrix
- `b::Array{Float, 1}`: right-hand side
- `c::Array{Float, 1}`: objective
- `tol`: numerical tolerance
- `verbose`: 0 means no output, 1 display logs at each iteration

"""
function solve(
    A,
    b,
    c,
    tol=10.0^-8,
    verbose=0
)
    # extract dimension
    m, n = size(A)
    @assert length(b) == m
    @assert length(c) == n

    # pre-allocate memory
    n_iter = 0
    mu = 0.0

    theta = zeros(n)
    b_norm = norm(b)
    c_norm = norm(c)

    # compute starting point
    x, y, s = startingPoint(A, b, c)

    rb = A * x - b   # primal residual
    rc = (y' * A + s' - c')'  # dual residual

    eps_p = (norm(rb)) / (1.0 + b_norm)  # relative primal feasibility
    eps_d = (norm(rc)) / (1.0 + b_norm)  # relative dual feasibility
    eps_g = abs(dot(c, x) - dot(b, y)) / (1.0 + dot(c, x))  # relative gap

    # X, Y and S are for analysis purposes
    # can be deleted to increase speed (and memory)
    X = [copy(x)]
    Y = [copy(y)]
    S = [copy(s)]

    N_ITER_MAX = 100  # maximum number of iterations

    println(" Itn      Primal Obj        Dual Obj    Prim Inf    Dual Inf\n")

    # main loop
    while (
        (
            (eps_p > tol)  # primal feasibility
            || (eps_d > tol)  # dual feasibility
            || (eps_g > tol)  # optimality gap
        )
        && (n_iter < N_ITER_MAX)
    )

        # book-keeping
        n_iter += 1
        @assert minimum(x) > 0
        @assert minimum(s) > 0

        mu = dot(x, s) / n  # centrality parameter
        theta = (x ./ s)
        
        # I. Form normal equations
        if n_iter == 1
            F = factorNormal(A, theta)
        else
            factorNormal!(F, A, theta)  # re-use factor if available
        end

        # II.1 Compute predictor step
        dx_aff, dy_aff, ds_aff = solveNormal(A, x, s, theta, -rc, -rb, -x.*s, F)

        # compute step size
        if minimum(dx_aff) >= 0.0
            a_prim = 1.0
        else
            a_prim = clamp(minimum(-(x ./dx_aff)[dx_aff .< 0]), 0.0, 1.0)
        end

        if minimum(ds_aff) >= 0.0
            a_dual = 1.0
        else
            a_dual = clamp(minimum(-(s ./ds_aff)[ds_aff .< 0]), 0.0, 1.0)
        end

        # update centrality parameter
        mu_aff = dot(
            x + a_prim * dx_aff,
            s + a_dual * ds_aff
        ) / n

        # II.2 Compute corrector step
        sigma = clamp((mu_aff / mu)^3, 10.0^-12, 1.0 - 10.0^-12)  # clamped for numerical stability

        ## corrector direction.
        dx_cc, dy_cc, ds_cc = solveNormal(
            A, x, s, theta,
            zeros(n), zeros(m),
            sigma*mu*ones(n) - dx_aff .* ds_aff,
            F
        )

        # compute step size
        dx = dx_aff + dx_cc
        dy = dy_aff + dy_cc
        ds = ds_aff + ds_cc

        # compute step size
        if minimum(dx) >= 0.0
            a_prim = 1.0
        else
            a_prim = clamp(minimum(-(x ./dx)[dx .< 0]), 0.0, 1.0)
        end

        if minimum(ds) >= 0.0
            a_dual = 1.0
        else
            a_dual = clamp(minimum(-(s ./ds)[ds .< 0]), 0.0, 1.0)
        end

        a_prim = min(0.99995 * a_prim, 1.0)
        a_dual = min(0.99995 * a_dual, 1.0)

        # III. Take step
        x += a_prim * dx
        y += a_dual * dy
        s += a_dual * ds

        push!(X, copy(x))
        push!(Y, copy(y))
        push!(S, copy(s))

        # compute residuals and stopping criterion
        rc = A' * y + s - c
        rb = A * x - b
        eps_p = (norm(rb)) / (1. + b_norm)
        eps_d = (norm(rc)) / (1. + c_norm)
        eps_g = abs(dot(c, x) - dot(b, y)) / (1. + abs(dot(c, x)))

        if verbose == 1
            print(@sprintf("%4.0f", n_iter))  # iteration count
            print(@sprintf("%+16.7e", dot(c, x)))  # primal objective
            print(@sprintf("%+16.7e", dot(b, y)))  # dual objective
            print(@sprintf("%+12.3e", norm(rb)))  # primal infeas
            print(@sprintf("%+12.3e", norm(rc)))  # dual infeas
            print("\n")
        end

    end

    # return solution
    return x, y, s, X, Y, S
end


"""
factorNormal(
    [F, ]
    A, theta
)

Form and factor normal equations' matrix.

# Arguments
- [F, ]: Existing factorization
- A: Matrix
- theta: vector
"""
function factorNormal(A, theta)
    # compute normal equations' matrix
    Phi = Symmetric((A * spdiagm(theta)) * A')

    # Cholesky factorization
    F = ldltfact(Phi)
    return F
end

function factorNormal!(F, A, theta)

    # compute normal equations' matrix
    Phi = Symmetric((A * spdiagm(theta)) * A')

    # Cholesky factorization
    ldltfact!(F, Phi)
    
    return nothing
end


"""
solveNormal(
    A, x, s, theta, xi_d, xi_p, xi_mu, F
)

Solve normal equations using existing factorization

# Arguments
- `A`:
- `x`:
- `s`:
- `theta`:
- `xi_d`:
- `xi_p`:
- `xi_mu`:
- `F`: Pre-computed factorization
"""
function solveNormal(
    A,
    x,
    s,
    theta,
    xi_d,
    xi_p,
    xi_mu,
    F
)

    # solve using the existing factorization
    dy = F \ (xi_p + A * (- xi_mu ./ s + theta .* xi_d))

    # back substitutions
    dx = theta .* (xi_mu ./ x - xi_d + (dy' * A)')
    ds = (xi_mu - s .* dx) ./ x

    return dx, dy, ds
end


"""
startingPoint(A, b, c)

Compute a (supposedly) good starting point for the IPM algorithm.

This is based on [Mehrotra 1992, On the implementation of a Primal-Dual
    Interior-point method], section 7.
"""
function startingPoint(A, b, c)

    # compute and factor A*A'
    Phi = A * A'
    F = cholfact(Symmetric(Phi))

    # compute initial feasible point
    y = F \ (A * c)
    s = (c' - y'*A)'
    x = F \ b
    x = (x' * A)'

    # compute deltas
    dx = max(-1.5*minimum(x), 0.0)
    ds = max(-1.5*minimum(s), 0.0)
    dx, ds = (
        dx + 0.5 * dot(x+dx, s+ds) / sum(s + ds),
        ds + 0.5 * dot(x+dx, s+ds) / sum(s + dx)
    )

    # compute starting point
    x += dx
    s += ds

    return x, y, s
end