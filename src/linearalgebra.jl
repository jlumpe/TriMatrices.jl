"""
	$(FUNCTIONNAME)(m::TriMatrix)

Wrap a [`TriMatrix`](@ref) in the appropriate `LinearAlgebra` view type for
optimized performance.
"""
function wraptri end

wraptri(m::TriUpperMatrix) = UpperTriangular(m)
wraptri(m::TriLowerMatrix) = LowerTriangular(m)
wraptri(m::TriSymmetricMatrix) = Symmetric(m)


Base.transpose(m::TriULMatrix) = TriMatrix(transpose_layout(TriLayout(m)), m.n, m.data; diag=m.diag)
Base.transpose(m::TriSymmetricMatrix) = m


LinearAlgebra.istriu(::TriUpperMatrix) = true
LinearAlgebra.istril(::TriLowerMatrix) = true
LinearAlgebra.issymmetric(::TriSymmetricMatrix) = true
