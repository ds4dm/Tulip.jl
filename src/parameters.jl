"""
    Parameters{T}

"""
Base.@kwdef mutable struct Parameters{T}
    # Model-wise parameters
    Threads::Int = 1
    OutputLevel::Int = 0

    # Linear algebra (MatrixOptions)
    MatrixFactory::Factory{<:AbstractMatrix} = Factory(SparseMatrixCSC)

    # Presolve
    Presolve::PresolveOptions{T} = PresolveOptions{T}()

    # IPM
    IPM::IPMOptions{T} = IPMOptions{T}()
    
    # KKT
    KKT::KKTOptions{T} = KKTOptions{T}()
end