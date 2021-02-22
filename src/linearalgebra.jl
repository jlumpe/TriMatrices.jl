"""
	$(FUNCTIONNAME)(m::TriMatrix)

Wrap a [`TriMatrix`](@ref) in the appropriate `LinearAlgebra` view type for
optimized performance.
"""
function wraptri end

wraptri(m::TriMatrix{<:TriLower}) = LowerTriangular(m)
wraptri(m::TriMatrix{<:TriUpper}) = UpperTriangular(m)
wraptri(m::TriMatrix{<:TriSymmetric}) = Symmetric(m)


Base.transpose(m::TriMatrix{TriUpper, L}) where L = TriMatrix{TriLower, L}(m.n, m.data, m.diag)
Base.transpose(m::TriMatrix{TriLower, L}) where L = TriMatrix{TriUpper, L}(m.n, m.data, m.diag)
Base.transpose(m::TriMatrix{TriSymmetric}) = m


LinearAlgebra.istriu(::TriMatrix{<:TriUpper}) = true
LinearAlgebra.istril(::TriMatrix{<:TriLower}) = true
LinearAlgebra.issymmetric(::TriMatrix{<:TriSymmetric}) = true
