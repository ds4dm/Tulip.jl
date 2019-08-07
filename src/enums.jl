"""
    BoundType

"""
@enum BoundType begin
    TLP_BND_UP  # Upper-bounded
    TLP_BND_LO  # Lower-bounded
    TLP_BND_FX  # Fixed
    TLP_BND_FR  # Free
    TLP_BND_RG  # Range
end