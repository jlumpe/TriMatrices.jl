"""
$(TYPEDEF)

Abstract type which indicates which regions of a [`TriMatrix`](@ref) are stored in
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


const TriUL{D} = Union{TriUpper{D}, TriLower{D}}


"""
	$(FUNCTIONNAME)(::Type{<:TriLayout})::Bool
	$(FUNCTIONNAME)(layout::TriLayout)::Bool

Check if the given layout/layout type stores values along the diagonal.
"""
hasdiag(::Type{<:TriLayout{D}}) where D = D
hasdiag(::TriLayout{D}) where D = D


"""
	$(FUNCTIONNAME)(::Type{<:TriLayout}, n)::Integer
	$(FUNCTIONNAME)(layout::TriLayout, n)::Integer

Get the number of elements needed to store the data of a `TriMatrix`
with the given layout.
"""
nelems(::Type{L}, n::Integer) where {L<:TriLayout} = trinum(max(hasdiag(L) ? n : n - 1, 0))
nelems(layout::TriLayout, n::Integer) = nelems(typeof(layout), n)


"""
	$(FUNCTIONNAME)(::Type{<:TriLayout})
	$(FUNCTIONNAME)(layout::TriLayout)

Get the `TriLayout` of the transpose of a `TriMatrix` with the given layout.

This maps `TriUpper` and `TriLower` to each other, and `TriSymmetric` to itself.
"""
transpose_layout(::Type{TriUpper{D}}) where D = TriLower{D}
transpose_layout(::Type{TriLower{D}}) where D = TriUpper{D}
transpose_layout(::Type{TriSymmetric{D}}) where D = TriSymmetric{D}
transpose_layout(::Type{TriUpper}) = TriLower
transpose_layout(::Type{TriLower}) = TriUpper
transpose_layout(::Type{TriSymmetric}) = TriSymmetric
transpose_layout(::L) where {L <: TriLayout} = transpose_layout(L)()
