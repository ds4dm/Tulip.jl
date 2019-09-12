"""
    ProblemData{Tv<:Real, Ti<:Integer}

Pace-holder and interface for problem data.

Problem data is stored in canonical form
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

    # Objective sense
    # (true is minimize, false is maximize)
    obj_sense::ObjSense
    obj_const::Tv  # Constant objective offset

    # Coefficients of the constraint matrix
    coeffs::Dict{Tuple{VarId, ConstrId}, Tv}
    var2con::Dict{VarId, OrderedSet{ConstrId}}
    con2var::Dict{ConstrId, OrderedSet{VarId}}

    # Variables
    vars::OrderedDict{VarId, Variable{Tv}}
    name2var::Dict{String, VarId}  # Matches variable name to index

    # Constraints
    constrs::OrderedDict{ConstrId, LinearConstraint{Tv}}
    name2con::Dict{String, ConstrId}  # Matches constraint name to index

    # Only allow empty problems to be instantiated for now
    function ProblemData{Tv}() where {Tv<:Real}
        return new{Tv}(
            0, 0, TLP_MIN, zero(Tv),
            Dict{Tuple{VarId, ConstrId}, Tv}(),
            Dict{VarId, OrderedSet{ConstrId}}(),
            Dict{ConstrId, OrderedSet{VarId}}(),
            OrderedDict{VarId, Variable{Tv}}(),
            Dict{String, VarId}(),
            OrderedDict{ConstrId, LinearConstraint{Tv}}(),
            Dict{String, ConstrId}()
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
    add_variable(pb::ProblemData{Tv}, v::Variable{Tv})

Add variable `v` to problem `pb`. Raises an error if `v` is already in the model. 
"""
function add_variable!(pb::ProblemData{Tv}, v::Variable{Tv}) where {Tv<:Real}
    !haskey(pb.vars, v.id) || error("Variable $(v.id.uuid) already exists.")
    name = v.dat.name
    # Ignore name if empty
    (name != "") || !haskey(pb.name2var, name) || error("Variable name $(name) already exists.")

    pb.vars[v.id] = v
    pb.var2con[v.id] = OrderedSet{ConstrId}()
    (name != "") && (pb.name2var[name] = v.id)
    return nothing
end


"""
    delete_variable!(pb::ProblemData{Tv}, vid::VarId)

Delete variable `vid` from the problem.
"""
function delete_variable!(pb::ProblemData{Tv}, vid::VarId) where{Tv<:Real}
    var = pb.vars[vid]  # to get variable name

    # Delete this column's coefficients
    cons = pb.var2con[vid]
    for cid in cons
        delete!(pb.coeffs, (vid, cid))
    end
    # Delete links in corresponding rows
    for cid in cons
        cols = pb.con2var[cid]
        delete!(cols, vid)
    end

    # Delete variable from var2con
    delete!(pb.var2con, vid)

    # Delete variable from name2var
    # Nothing to delete if name was empty
    if var.dat.name != ""
        delete!(pb.name2var, var.dat.name)
    end

    # Finally, delete variable from list of variables
    delete!(pb.vars, vid)

    return nothing
end


"""
    add_constraint!(pb::ProblemData, c::LinearConstraint{Tv})

Create a new linear constraint, add it to the model, and return its ID.
"""
function add_constraint!(pb::ProblemData{Tv}, c::LinearConstraint{Tv}) where{Tv<:Real}
    !haskey(pb.constrs, c.id) || error("Constraint $(c.id.uuid) already exists.")
    name = c.dat.name
    # Ignore name if empty
    (name == "") || !haskey(pb.name2var, name) || error("Variable name $(name) already exists.")

    pb.constrs[c.id] = c
    pb.con2var[c.id] = OrderedSet{VarId}()
    (name != "") && (pb.name2con[name] = c.id)
    return nothing
end


"""
    delete_constraint!(pb::ProblemData{Tv}, cid::ConstrId)

Delete constraint `cid` from the problem.
"""
function delete_constraint!(pb::ProblemData{Tv}, cid::ConstrId) where{Tv<:Real}
    con = pb.constrs[cid]  # to get constraint name

    # Delete this row's coefficients
    vars = pb.con2var[cid]
    for vid in vars
        delete!(pb.coeffs, (vid, cid))
    end
    # Delete links in corresponding columns
    for vid in vars
        rows = pb.var2con[vid]
        delete!(rows, cid)
    end

    # Delete constraint from con2var
    delete!(pb.con2var, cid)

    # Delete variable from name2con
    # Nothing to delete if name was empty
    if con.dat.name != ""
        delete!(pb.name2con, con.dat.name)
    end

    # Finally, delete constraint from list of constraints
    delete!(pb.constrs, cid)

    return nothing
end


"""
    set_coeff!(pb::ProblemData{Tv}, vid::VarId, cid::ConstrId, val::Tv)

"""
function set_coeff!(pb::ProblemData{Tv}, vid::VarId, cid::ConstrId, val::Real) where{Tv<:Real}
    
    if !haskey(pb.coeffs, (vid, cid))
        # Entry does not exist
        iszero(val) && return nothing  # Ignore zero values

        # Set coefficient
        pb.coeffs[vid, cid] = val
        push!(pb.var2con[vid], cid)
        push!(pb.con2var[cid], vid)
    else
        # Entry exists. Need to check if coeff is zero or not
        if iszero(val)
            # Delete coefficient from problem
            delete!(pb.coeffs, (vid, cid))
            delete!(pb.var2con[vid], cid)
            delete!(pb.con2var[cid], vid)
        else
            # Just replace the value
            pb.coeffs[vid, cid] = val
        end
    end

    return nothing
end