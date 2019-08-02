"""
    ProblemData{Tv<:Real, Ti<:Integer}

Problem data is stored in the form
```math
\\begin{align}
    \\min_{x} \\ \\ \\ & c^{T} x \\\\
    s.t. \\ \\ \\
    & l_c \\leq A x \\leq u_c \\\\
    & l_x \\leq x \\leq u_x
\\end{align}
```

and will be converted to standard form
    min     c'x
    s.t.    A x = b
            0 <= x <= u
before the optimizer is called.

Each variable and constraint has a unique index.
"""
mutable struct ProblemData{Tv<:Real, Ti<:Integer}

    ncon::Ti  # Number of constraints
    nvar::Ti  # Number of variables

    # Coefficients of the constraint matrix
    coeffs::Dict{Tuple{Ti, Ti}, Tv}
    var2con::Dict{Ti, OrderedSet{Ti}}
    con2var::Dict{Ti, OrderedSet{Ti}}

    # Variables
    vars::OrderedDict{Ti, Variable{Tv, Ti}}

    # Constraints
    constrs::OrderedDict{Ti, AbstractConstraint{Tv, Ti}}

    # Only allow empty problems to be instantiated for now
    function ProblemData{Tv, Ti}() where {Tv<:Real, Ti<:Integer}
        return new{Tv, Ti}(
            zero(Ti), zero(Ti),
            Dict{Tuple{Ti, Ti}, Tv}(),
            Dict{Ti, OrderedSet{Ti}}(), Dict{Ti, OrderedSet{Ti}}(),
            OrderedDict{Ti, Variable{Tv, Ti}}(),
            OrderedDict{Ti, AbstractConstraint{Tv, Ti}}()
        )
    end
end


"""
    get_num_var(pb)

Return the number of variables in the problem.
"""
get_num_var(pb::ProblemData) = length(pb.vars)


"""
    get_num_constr(pb)

Return the number of constraints in the problem.
"""
get_num_constr(pb::ProblemData) = length(pb.constrs)


# TODO: replace pb.nvar += 1 by an increment function that checks if typemax
# is reached, and raises an error if it is.
"""
    add_variable(pb::ProblemData)

Create a new variable in the model and return the corresponding ID.
"""
function add_variable!(pb::ProblemData{Tv, Ti}) where {Tv, Ti}
    uuid = pb.nvar + 1
    pb.nvar += 1

    v = Variable{Tv, Ti}(uuid)
    pb.vars[uuid] = v
    return uuid
end


"""
    add_linear_constraint(pb::ProblemData)

Create a new linear constraint, add it to the model, and return its ID.
"""
function add_linear_constraint!(pb::ProblemData{Tv, Ti}) where {Tv, Ti}
    uuid = pb.ncon + 1
    pb.ncon += 1

    c = LinearConstraint{Tv, Ti}(uuid)
    pb.constrs[uuid] = c
    return uuid
end