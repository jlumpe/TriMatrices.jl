using SafeTestsets


@safetestset "Testing" begin include("testing.jl") end
@safetestset "Layouts" begin include("layouts.jl") end
@safetestset "Indexing" begin include("indexing.jl") end
@safetestset "TriMatrix" begin include("trimatrix.jl") end
@safetestset "LinearAlgebra" begin include("linearalgebra.jl") end
