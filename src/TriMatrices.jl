module TriMatrices

using DocStringExtensions

export TriUpper, TriLower, TriSymmetric


include("layout.jl")

include("Indexing/Indexing.jl")
using .Indexing

include("Testing.jl")


end # module
