"""
    Residuals{T, Tv}

Data structure for IPM residual vectors.
"""
mutable struct Residuals{T, Tv}
    # Primal residuals
    rp::Tv  # rp = τ*b - A*x
    rl::Tv  # rl = τ*l - (x - xl)
    ru::Tv  # ru = τ*u - (x + xu)

    # Dual residuals
    rd::Tv  # rd = τ*c - (A'y + zl - zu)
    rg::T  # rg = c'x - (b'y + l'zl - u'zu) + κ

    # Residuals' norms
    rp_nrm::T  # |rp|
    rl_nrm::T  # |rl|
    ru_nrm::T  # |ru|
    rd_nrm::T  # |rd|
    rg_nrm::T  # |rg|
end