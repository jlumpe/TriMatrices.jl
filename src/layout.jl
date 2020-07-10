"""
$(TYPEDEF)

Abstract type which indicates how a square matrix's values are layed out in
memory. Type parameter is a `Bool` indicating whether the matrix diagonal is
stored.
"""
abstract type TriLayout{D} end

# Check the type parameter of a TriLayout and throw an error if it isn't a Bool
_check_layout_param(D) = D isa Bool || throw(ArgumentError("Type parameter of TriLayout must be a Bool"))

# Define a new TriLayout subtype
macro trilayout(T::Symbol)
	return quote
		struct $T{D} <: TriLayout{D}
			$T{D}() where D = (_check_layout_param(D); new{D}())
			$T() = $T{true}()
		end
	end |> esc
end

"""
$(TYPEDEF)

The upper triangle of the matrix is stored, values beneath the diagonal are
zero.
"""
@trilayout TriUpper

"""
$(TYPEDEF)

The lower triangle of the matrix is stored, values above the diagonal are
zero.
"""
@trilayout TriLower

"""
$(TYPEDEF)

Matrix is symmetric across the diagonal, one number is stored for each pair of
non-diagonal entries.
"""
@trilayout TriSymmetric


"""
	$(FUNCTIONNAME)(layout::TriLayout)

Check if the given layout stores values along the diagonal.
"""
hasdiag(::Type{<:TriLayout{D}}) where D = D
hasdiag(::TriLayout{D}) where D = D


"""
	$(FUNCTIONNAME)(layout::TriLayout, n::Integer)

Get the number of elements needed to store the data of an `n` by `n` matrix
with the given layout.
"""
nelems(layout, n::Integer) = trinum(hasdiag(layout) ? n : n - 1)
