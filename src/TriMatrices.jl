module TriMatrices

using DocStringExtensions

export TriUpper, TriLower, TriSymmetric


include("math.jl")
include("layout.jl")
include("indexing.jl")
include("testing.jl")


end # module
