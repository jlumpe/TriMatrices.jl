module Indexing

using DocStringExtensions
using TriMatrices: TriLayout, TriLower, TriUpper, TriSymmetric,
                   nelems, hasdiag


include("math.jl")
include("conversion.jl")
include("tri_indices.jl")

end  # module
