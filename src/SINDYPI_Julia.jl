module SINDYPI_Julia

include("SINDYPI_main.jl")
include("Basis_generation.jl")
include("SparsityAlgorithm.jl")
include("Table_generation.jl")

export  SINDYPI
export Basis_SINDY_PI
export Sparse_matrix_SINDYPI
export Table_SINDY_PI

end # module
