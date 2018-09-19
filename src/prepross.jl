function prepross!(m::Model)

    # Right-hand side
    m.b = zeros(m.num_constr)

    # Process constraints
    for i in Base.OneTo(m.num_constr)

        l = m.constr_lb[i]
        u = m.constr_ub[i]

        # 
        if l == -Inf && u == Inf
            # Constraint is always satisfied, raise error
            error("Constraint $i is always true")
        
        elseif l > u || l == Inf || u == -Inf
            # Constraint is infeasible
            error("Constraint $i is infeasible")

        elseif l == -Inf
            # Lesser-than inequality
            # a'x ≦ u
            # Add slack variable
            addvar!(m, SparseVector(m.num_constr, [i], [1.0]), 0.0 , Inf, 0.0)
            m.b[i] = u

        elseif u == Inf
            # Greater-than inequality
            # Add slack variable
            addvar!(m, SparseVector(m.num_constr, [i], [-1.0]), 0.0 , Inf, 0.0)
            m.b[i] = l

        elseif -Inf < l < u < Inf
            # Range constraint
            # Add slack variable with lower and upper bound
            addvar!(m, SparseVector(m.num_constr, [i], [-1.0]), l , u, 0.0)
            m.b[i] = 0.0
        else
            # Equality constraint
            # No transformation needed
            m.b[i] = l
        end
        
    end

    # Process variables
    m.c = zeros(m.num_var)
    m.uind = Vector{Int}(undef, 0)
    m.uval = zeros(0)

    for i in Base.OneTo(m.num_var)
        
        l = m.var_lb[i]
        u = m.var_ub[i]

        if l > u || l == Inf || u == -Inf
            # Bounds constraints cannot be satisfied
            error("Bounds constraints on variable $i are infeasible.")
        
        elseif l == -Inf && u == Inf
            # Free variable
            error("Free variables are not supported.")

        elseif l == -Inf && u < Inf
            # Upper bounded variable, no lower bound
            m.c[i] = -m.obj[i]
            m.A[:, i] .*= -1

            if u == 0.0
                # New variable is non-negative
            else
                # Offset right-hand side
                axpy!(u, m.A[:, i], m.b)
            end

        elseif l > -Inf && u == Inf
            m.c[i] = m.obj[i]

            # Lower-bounded variable, no upper bound
            if l != 0.0
                # Offset right-hand side
                axpy!(-l, m.A[:, i], m.b)
            end
        else 
            m.c[i] = m.obj[i]
            # Bounded variable
            push!(m.uind, i)
            push!(m.uval, u-l)

            # Offset right-hand side
            if l != 0.0
                # Offset right-hand side
                axpy!(-l, m.A[:, i], m.b)
            end
        end
        
    end
    m.n_var_ub = length(m.uind)

end