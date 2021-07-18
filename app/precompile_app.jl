using TulipCL

const EXDIR = joinpath(dirname(pathof(TulipCL.Tulip)), "../examples/dat")

# Run all example problems
for finst in readdir(EXDIR)
    empty!(ARGS)
    append!(ARGS, ["--Threads", "1", "--TimeLimit", "10.0", "--Presolve", "1", "--Method", "HSD", joinpath(EXDIR, finst)])
    TulipCL.julia_main()
end

const NETLIB = TulipCL.Tulip.QPSReader.fetch_netlib()

for finst in readdir(NETLIB)[1:5], ipm in ["HSD", "MPC"]
    empty!(ARGS)
    append!(ARGS, ["--Threads", "1", "--TimeLimit", "10.0", "--Presolve", "1", "--Method", ipm, joinpath(NETLIB, finst)])
    TulipCL.julia_main()
end
