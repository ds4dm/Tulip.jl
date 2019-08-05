"""
    ProblemData{Tv<:Real, Ti<:Integer}

Pace-holder and interface for problem data.

Problem data is stored in the form
```math
\\begin{align}
    \\min_{x} \\ \\ \\ & c^{T} x \\\\
    s.t. \\ \\ \\
    & l_c \\leq A x \\leq u_c \\\\
    & l_x \\leq x \\leq u_x
\\end{align}
```
"""
mutable struct ProblemData{Tv<:Real}

    ncon::Int  # Number of constraints
    nvar::Int  # Number of variables

    # Coefficients of the constraint matrix
    coeffs::Dict{Tuple{VarId, ConstrId}, Tv}
    var2con::Dict{VarId, OrderedSet{ConstrId}}
    con2var::Dict{ConstrId, OrderedSet{VarId}}

    # Variables
    vars::OrderedDict{VarId, Variable{Tv}}

    # Constraints
    constrs::OrderedDict{ConstrId, AbstractConstraint{Tv}}

    # Only allow empty problems to be instantiated for now
    function ProblemData{Tv}() where {Tv<:Real}
        return new{Tv}(
            0, 0,
            Dict{Tuple{VarId, ConstrId}, Tv}(),
            Dict{VarId, OrderedSet{ConstrId}}(),
            con2var::Dict{ConstrId, OrderedSet{VarId}}(),
            OrderedDict{VarId, Variable{Tv}}(),
            OrderedDict{ConstrId, AbstractConstraint{Tv}}()
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
function add_variable!(pb::ProblemData{Tv}) where {Tv}
    uuid = pb.nvar + 1
    pb.nvar += 1

    v = Variable{Tv}(uuid)
    pb.vars[uuid] = v
    return uuid
end


"""
    add_linear_constraint(pb::ProblemData)

Create a new linear constraint, add it to the model, and return its ID.
"""
function add_linear_constraint!(pb::ProblemData{Tv}) where {Tv}
    uuid = pb.ncon + 1
    pb.ncon += 1

    c = LinearConstraint{Tv}(uuid)
    pb.constrs[uuid] = c
    return uuid
end