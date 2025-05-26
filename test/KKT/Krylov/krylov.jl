using Krylov

if !isdefined(Krylov, :CgWorkspace)
    const CgWorkspace = Krylov.CgSolver
    const MinresWorkspace = Krylov.MinresSolver
    const MinresQlpWorkspace = Krylov.MinresQlpSolver
    const SymmlqWorkspace = Krylov.SymmlqSolver
    const TricgWorkspace = Krylov.TricgSolver
    const TrimrWorkspace = Krylov.TrimrSolver
end

@testset "Krylov" begin
    include("spd.jl")
    include("sid.jl")
    include("sqd.jl")
end
