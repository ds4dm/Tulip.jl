using SparseArrays

"""
    StandardForm{Tv, Ti, Ta}

Problem data converted in standard form.
```math
\\begin{align}
    \\min_{x} \\ \\ \\ & c^{T} x \\\\
    s.t. \\ \\ \\
    & A x = b \\\\
    & x_{i} \\leq u \\\\
    & x \\geq 0
\\end{align}
```
"""
mutable struct StandardForm{Tv<:Real}
    A::AbstractMatrix{Tv}  # Constraint matrix
    b::Vector{Tv}          # Right-hand side
    c::Vector{Tv}          # Objective
    uind::Vector{Int}      # Indices of upper-bounded variables
    uval::Vector{Tv}       # Finite upper bounds on variables

    # TODO: add optimization sense
    # TODO: add row and column scalings
    # TODO: add starting points
    # TODO: add Cholesky factor (?)
end


"""
    convert_to_standard_form(pb::ProblemData{Tv})

Convert problem to standard form.
"""
function convert_to_standard_form(pb::ProblemData{Tv}) where {Tv<:Real}

    nvars = 0  # Number of variables   in standard form model
    ncons = 0  # Number of constraints in standard form model

    m0 = get_num_constr(pb)  # number of original variables
    n0 = get_num_var(pb)  # number of original constraints

    # Keep track of indices
    con2idx = Dict{ConstrId, Int}()
    var2idx = Dict{VarId, Int}()
    idx2var = Union{VarId, ConstrId}[]
    idx2con = ConstrId[]

    # Get right-hand side and allocate b
    b = zeros(Tv, m0)
    nslack = 0  # Number of additional slack variables
    for (cidx, con) in pb.constrs
        if con.dat.bt == TLP_FX
            # Equality constraint: keep as is
            ncons += 1
            con2idx[cidx] = ncons
            push!(idx2con, cidx)
            b[ncons] = con.dat.lb
        elseif con.dat.bt == TLP_LO
            # a'x >= b: add surplus
            ncons += 1
            nslack += 1
            con2idx[cidx] = ncons
            push!(idx2con, cidx)
            b[ncons] = con.dat.lb
        elseif con.dat.bt == TLP_UP
            # a'x <= b: add slack
            ncons += 1
            nslack += 1
            con2idx[cidx] = ncons
            push!(idx2con, cidx)
            b[ncons] = con.dat.ub
        elseif con.dat.bt == TLP_RG
            # Range constraint: add a bounded slack
            ncons += 1
            nslack += 1
            con2idx[cidx] = ncons
            push!(idx2con, cidx)
            b[ncons] = con.dat.lb
        else
            error("Unsupported bound type: $(con.dat.bt)")
        end
    end

    # Allocate objective and coefficients
    c = Tv[]
    aI = Int[]
    aJ = Int[]
    aV = Tv[]
    uind = Int[]
    uval = Tv[]

    for (vidx, var) in pb.vars
        # Record variable info

        if var.dat.bt == TLP_FX
            # Fixed variable
            # TODO: remove
            nvars += 1
            push!(c, var.dat.obj)
            var2idx[vidx] = nvars
            push!(idx2var, vidx)
            vind = nvars

            # Add upper bound
            push!(uind, vind)
            push!(uval, zero(Tv))

            for cidx in pb.var2con[vidx]
                # add coeffs to A
                v = pb.coeffs[vidx, cidx]
                cind = con2idx[cidx]

                push!(aI, cind)
                push!(aJ, vind)
                push!(aV, v)

                # Update riht-hand side
                b[cind] -= v * var.dat.lb
            end
            
        elseif var.dat.bt == TLP_UP
            # Upper-bounded variable: x <= u

            # Change to -x >= -u and shift right-hand side
            # Corresponding column and obj coeff are multiplied by -1
            nvars += 1
            push!(c, -var.dat.obj)
            var2idx[vidx] = nvars
            push!(idx2var, vidx)
    
            vind = nvars
            for cidx in pb.var2con[vidx]
                # add coeffs to A
                v = pb.coeffs[vidx, cidx]
                cind = con2idx[cidx]

                push!(aI, cind)
                push!(aJ, vind)
                push!(aV, -v)

                # Update right-hand side
                b[cind] += v * var.dat.ub
            end


        elseif var.dat.bt == TLP_LO
            # Lower-bounded variable x >= l

            # Shift right-hand side
            nvars += 1
            push!(c, var.dat.obj)
            var2idx[vidx] = nvars
            push!(idx2var, vidx)
            vind = nvars

            for cidx in pb.var2con[vidx]
                # add coeffs to A
                v = pb.coeffs[vidx, cidx]
                cind = con2idx[cidx]
                push!(aI, cind)
                push!(aJ, vind)
                push!(aV, v)

                # Update right-hand side
                b[cind] -= v * var.dat.lb
            end

        elseif var.dat.bt == TLP_FR
            # Free variable
            # Split into the difference of two variables

            # TODO: properly keep track of that
            nvars += 2
            push!(c,  var.dat.obj)  # positive part
            push!(c, -var.dat.obj)  # negative part
            var2idx[vidx] = nvars -1  # Just recall the first 
            push!(idx2var, vidx)
            push!(idx2var, vidx)

            vind = nvars - 1

            for cidx in pb.var2con[vidx]
                # add coeffs to A
                v = pb.coeffs[vidx, cidx]
                cind = con2idx[cidx]

                # Positive part
                push!(aI, cind)
                push!(aJ, vind)
                push!(aV, v)

                # Negative part
                push!(aI, cind)
                push!(aJ, vind + 1)
                push!(aV, -v)
            end

        elseif var.dat.bt == TLP_RG
            # Ranged variable
            # Shift right-hand side
            # update upper-bounds vectors
            nvars += 1
            push!(c, var.dat.obj)
            var2idx[vidx] = nvars
            push!(idx2var, vidx)
            vind = nvars

            for cidx in pb.var2con[vidx]
                # add coeffs to A
                v = pb.coeffs[vidx, cidx]
                cind = con2idx[cidx]

                push!(aI, cind)
                push!(aJ, vind)
                push!(aV, v)

                # Update right-hand side
                b[cind] -= v * var.dat.lb
            end

            # Update bounds
            push!(uind, vind)
            push!(uval, var.dat.ub - var.dat.lb)

        else
            error("Unsupported bound type $(var.dat.bt)")
        end
    end  # variables loop


    # Finally, add coefficients for the slacks
    # Go through the constraints once again
    # Here we need to update c, uind, uval and the coeffs of A 
    for (cidx, con) in pb.constrs

        cind = con2idx[cidx]

        if con.dat.bt == TLP_LO
            # a'x >= l
            # change to a'x - s = l
            nvars += 1
            vind = nvars
            push!(idx2var, cidx)

            # Update objective
            push!(c, zero(Tv))

            # Add coefficient
            push!(aI, cind)
            push!(aJ, vind)
            push!(aV, -oneunit(Tv))

        elseif con.dat.bt == TLP_UP
            # a'x <= u
            # change to a'x + s = u
            nvars += 1
            vind = nvars
            push!(idx2var, cidx)

            # Update objective
            push!(c, zero(Tv))

            # Add coefficient
            push!(aI, cind)
            push!(aJ, vind)
            push!(aV, oneunit(Tv))

        elseif con.dat.bt == TLP_RG
            # l <= a'x <= u
            # change to a'x - s = l, 0 <= s <= u-l
            nvars += 1
            vind = nvars
            push!(idx2var, cidx)

            # Update objective
            push!(c, zero(Tv))

            # Add coefficient
            push!(aI, cind)
            push!(aJ, vind)
            push!(aV, -oneunit(Tv))

            # Update upper bounds
            push!(uind, vind)
            push!(uval, con.dat.ub - con.dat.lb)

        end
    end

    # Done.

    
    return ncons, nvars, aI, aJ, aV, b, c, uind, uval, con2idx, var2idx, idx2con, idx2var
end