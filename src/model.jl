"""
    Model
    
Data structure for a model.

# Attributes
"""
mutable struct Model
    #=======================================================
        Model attributes
    =======================================================#

    env::TulipEnv                   # Environment
    sln_status::SolutionStatus      # Solution status
    
    time_total::Float64     # Total time spent by the optimization
    num_bar_iter::Int       # Total number of Interior-Point iterations
    primal_bound::Float64   # Best known primal bound
    dual_bound::Float64     # Best known dual bound


    #=======================================================
        Problem data
    =======================================================#

    num_var::Int              # Number of variables (total)
    n_var_ub::Int             # Number of upper-bounded variables
    num_constr::Int           # Number of constraints

    A::AbstractMatrix{Float64}      # Constraints matrix
    obj::AbstractVector{Float64}    # Objective coefficients
    var_lb::Vector{Float64}         # Lower bounds on `x` variables
    var_ub::Vector{Float64}         # Upper bounds on `x` variables
    constr_lb::Vector{Float64}      # Lower bounds on linear constraints
    constr_ub::Vector{Float64}      # Upper bounds on linear constraints


    #=======================================================
        Working memory
    =======================================================#
    b::AbstractVector{Float64}      # Right-hand side of equality constraints
    c::AbstractVector{Float64}      # Objective vector
    uind::AbstractVector{Int}       # Indices of upper-bounded variables
    uval::AbstractVector{Float64}   # Values of upper bounds

    x::Vector{Float64}      # Vector of primal variables (original variables)
    w::Vector{Float64}      # Vector of primal variables (upper bound slack)
    y::Vector{Float64}      # Vector of dual variables (equality constraints)
    s::Vector{Float64}      # Vector of reduced costs of `x`
    z::Vector{Float64}      # Vector of reduced costs of `w`
    t::Base.RefValue{Float64}    # Artificial homogeneous primal variable
    k::Base.RefValue{Float64}    # Artificial homogeneous dual variable
    μ::Base.RefValue{Float64}    # Current barrier parameter

    rp::Vector{Float64}     # Vector of primal residuals `Ax - b`
    rd::Vector{Float64}     # Vector of dual residuals `A'y + s - c`
    ru::Vector{Float64}     # Vector of primal residuals 'x + w - u`
    rg::Base.RefValue{Float64}   # Residual for optimality gap
    rxs::Vector{Float64}    # Right-hand side for complimentarity product `x*s`
    rwz::Vector{Float64}    # Right-hand side for complimentarity product `w*z`
    rtk::Base.RefValue{Float64}  # Right-hand side for complimentarity product `t*k`


    #=======================================================
        Model constructor
    =======================================================#  
    function Model(
        env::TulipEnv
    )
        m = new()
        m.env = env
        
        # Initialize data memory
        m.num_var = 0
        m.num_constr = 0
        m.A = spzeros(0, 0)
        m.obj = Vector{Float64}(undef, 0)
        m.var_lb = Vector{Float64}(undef, 0)
        m.var_ub = Vector{Float64}(undef, 0)
        m.constr_lb = Vector{Float64}(undef, 0)
        m.constr_ub = Vector{Float64}(undef, 0)

        # Initialize working memory
        m.x = Vector{Float64}(undef, 0)
        m.w = Vector{Float64}(undef, 0)
        m.y = Vector{Float64}(undef, 0)
        m.s = Vector{Float64}(undef, 0)
        m.z = Vector{Float64}(undef, 0)
        m.t = Ref(1.0)
        m.k = Ref(1.0)
        m.μ = Ref(1.0)

        m.rp = Vector{Float64}(undef, 0)
        m.ru = Vector{Float64}(undef, 0)
        m.rd = Vector{Float64}(undef, 0)
        m.rg = Ref(Inf)
        m.rxs = Vector{Float64}(undef, 0)
        m.rwz = Vector{Float64}(undef, 0)
        m.rtk = Ref(Inf)

        # Book-keeping stuff
        m.sln_status = Sln_Unknown

        return m
    end
end

include("model_api.jl")