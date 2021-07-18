module TulipCL

using Printf

import Tulip
const TLP = Tulip

using ArgParse

function julia_main()::Cint
    try
        tulip_cl()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function parse_commandline(cl_args)

    s = ArgParseSettings()

    @add_arg_table! s begin
        "--TimeLimit"
            help = "Time limit, in seconds."
            arg_type = Float64
            default = Inf
        "--IterationsLimit"
            help = "Maximum number of iterations"
            arg_type = Int
            default = 500
        "--Threads"
            help = "Maximum number of threads."
            arg_type = Int
            default = 1
        "--Presolve"
            help = "Presolve level"
            arg_type = Int
            default = 1
        "--Method"
            help = "Interior-point method (HSD or MPC)"
            arg_type = String
            default = "HSD"
        "finst"
            help = "Name of instance file. Only Free MPS format is supported."
            required = true
    end

    return parse_args(cl_args, s)
end

function tulip_cl()
    parsed_args = parse_commandline(ARGS)

    # Read model and solve
    finst::String = parsed_args["finst"]

    m = TLP.Model{Float64}()
    t = @elapsed TLP.load_problem!(m, finst)

    println("Julia version: ", VERSION)
    println("Tulip version: ", Tulip.version())
    println("Problem file : ", finst)
    @printf("Reading time : %.2fs\n\n", t)

    # Set parameters
    m.params.OutputLevel = 1
    m.params.IPM.TimeLimit = parsed_args["TimeLimit"]
    m.params.Threads = parsed_args["Threads"]
    m.params.Presolve.Level = parsed_args["Presolve"]
    m.params.IPM.IterationsLimit = parsed_args["IterationsLimit"]

    if parsed_args["Method"] == "HSD"
        m.params.IPM.Factory = Tulip.Factory(Tulip.HSD)
    elseif parsed_args["Method"] == "MPC"
        m.params.IPM.Factory = Tulip.Factory(Tulip.MPC)
    else
        error("Invalid value for Method: $(parsed_args["Method"]) (must be HSD or MPC)")
    end

    TLP.optimize!(m)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    tulip_cl()
end

end # module
