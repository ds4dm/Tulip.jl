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

    constr_cnt::Int  # Counter for constraints
    var_cnt::Int     # Counter for variables

    # Coefficients of the constraint matrix
    coeffs::Dict{Tuple{VarId, ConstrId}, Tv}
    var2con::Dict{VarId, OrderedSet{ConstrId}}
    con2var::Dict{ConstrId, OrderedSet{VarId}}

    # Variables
    vars::OrderedDict{VarId, Variable{Tv}}

    # Constraints
    constrs::OrderedDict{ConstrId, LinearConstraint{Tv}}

    # Only allow empty problems to be instantiated for now
    function ProblemData{Tv}() where {Tv<:Real}
        return new{Tv}(
            0, 0,
            Dict{Tuple{VarId, ConstrId}, Tv}(),
            Dict{VarId, OrderedSet{ConstrId}}(),
            Dict{ConstrId, OrderedSet{VarId}}(),
            OrderedDict{VarId, Variable{Tv}}(),
            OrderedDict{ConstrId, LinearConstraint{Tv}}()
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
    add_variable(pb::ProblemData{Tv}, v::Variable{Tv})

Add variable `v` to problem `pb`. Raises an error if `v` is already in the model. 
"""
function add_variable!(pb::ProblemData{Tv}, v::Variable{Tv}) where {Tv<:Real}
    !haskey(pb.vars, v.id) || error("Variable $(v.id.uuid) already exists.")

    pb.vars[v.id] = v
    pb.var2con[v.id] = OrderedSet{ConstrId}()
    return nothing
end


"""
    new_variable_index!(pb::ProblemData)

Returns a new variable index and update internal counter.
"""
function new_variable_index!(pb::ProblemData)
    idx = VarId(pb.var_cnt + 1)
    pb.var_cnt += 1
    return idx
end



"""
    new_constraint_index!(pb::ProblemData)

Returns a new constraint index and update internal counter.
"""
function new_constraint_index!(pb::ProblemData)
    idx = ConstrId(pb.constr_cnt + 1)
    pb.constr_cnt += 1
    return idx
end


"""
    add_constraint!(pb::ProblemData, c::LinearConstraint{Tv})

Create a new linear constraint, add it to the model, and return its ID.
"""
function add_constraint!(pb::ProblemData{Tv}, c::LinearConstraint{Tv}) where{Tv<:Real}
    !haskey(pb.constrs, c.id) || error("Constraint $(c.id.uuid) already exists.")
    
    pb.constrs[c.id] = c
    pb.con2var[c.id] = OrderedSet{VarId}()
    return nothing
end