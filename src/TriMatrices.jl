module TriMatrices

using LinearAlgebra
using DocStringExtensions

export TriUpper, TriLower, TriSymmetric, TriMatrix


include("layout.jl")

include("Indexing/Indexing.jl")
using .Indexing

include("trimatrix.jl")
include("linearalgebra.jl")
include("Testing.jl")


end # module
