@testset "Dense" begin
    
# SparseCSC matrix
A = [
    1.0 0.0 1.0 0.0;
    0.0 1.0 0.0 1.0
]
m, n = size(A)
c = [1.1, 1.2, 1.3, 1.4]
b = [1.1, 0.9]
u = sparse([0.5, 0.4, 0.7, 1.1])
uind = u.nzind
uval = u.nzval
p = nnz(u)

test_linalg(
    A, b, c, uind, uval,
    zeros(n), zeros(p), zeros(m), zeros(n), zeros(p)
)

end  # testset