using SparseArrays

function presolve_test(Tv::Type)
    m = Tulip.Model{Tv}()

    # Build the following model
    #=
        min     x1 + x2 + x3
        s.t.    x1 + x2      = 1
                     x2 + x3 = 1
                x1,  x2,  x3 ⩾ 0
    =#
    pb = Tulip.ProblemData{Tv}()

    A = sparse(Tv.([
        [1 1 0];
        [0 1 1]
    ]))
    b = ones(Tv, m)
    c = ones(Tv, n)

    Tulip.load_problem!(pb, "test",
        true, c, zero(Tv),
        A, b, b, zeros(Tv, n), fill(Tv(Inf), n),
        ["c1", "c2"], ["x1", "x2", "x3"]
    )

    ps = Tulip.PresolveData(pb)

    return
end

include("./empty_column.jl")
include("./empty_row.jl")
include("./fixed_variable.jl")
include("./row_singleton.jl")

# for Tv in TvTYPES
#     @testset "$Tv" begin presolve_test(Tv) end
# end