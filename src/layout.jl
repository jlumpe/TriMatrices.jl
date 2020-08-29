"""
$(TYPEDEF)

Abstract type which indicates how a [`TriMatrix`](@ref)'s values are layed out in
memory. Type parameter is a `Bool` indicating whether the matrix diagonal is
stored.
"""
abstract type TriLayout{D} end

# Check the type parameter of a TriLayout and throw an error if it isn't a Bool
_check_layout_param(D) = D isa Bool || throw(ArgumentError("Type parameter of TriLayout must be a Bool"))


for T in [:TriUpper, :TriLower, :TriSymmetric]
	@eval struct $T{D} <: TriLayout{D}
		$T{D}() where D = (_check_layout_param(D); new{D}())
		$T() = $T{true}()
	end
end

@doc """
$(TYPEDEF)

The upper triangle of the matrix is stored, values beneath the diagonal are
zero.
"""
TriUpper

@doc """
$(TYPEDEF)

The lower triangle of the matrix is stored, values above the diagonal are
zero.
"""
TriLower

@doc """
$(TYPEDEF)

Matrix is symmetric across the diagonal, one value is stored for each pair of
non-diagonal entries.
"""
TriSymmetric


"""
	$(FUNCTIONNAME)(::Type{<:TriLayout})
	$(FUNCTIONNAME)(layout::TriLayout)

Check if the given layout/layout type stores values along the diagonal.
"""
hasdiag(::Type{<:TriLayout{D}}) where D = D
hasdiag(::TriLayout{D}) where D = D


"""
$(TYPEDSIGNATURES)

Get the number of elements needed to store the data of a `TriMatrix`
with the given layout.
"""
nelems(layout, n::Integer) = trinum(hasdiag(layout) ? n : n - 1)
