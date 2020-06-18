export tri_indices


"""
$(TYPEDEF)

Iterator over Cartesian indices of the stored region of a [`TriMatrix`](@ref), in
the same order as the corresponding in its data array.
"""
struct TriIndexIterator{L<:TriLayout}
	n::Int
end

Base.length(ti::TriIndexIterator{L}) where L = nelems(L, ti.n)
TriLayout(::TriIndexIterator{L}) where L = L()


Base.iterate(e::TriIndexIterator) = iterate(e, (1, 1))

function Base.iterate(e::TriIndexIterator{L}, state::Tuple{Int, Int}) where L
	i, j = state

	# Reached end of column
	if j > i
		i += 1
		j = 1
	end

	# Finished last row
	i > (hasdiag(L) ? e.n : e.n - 1) && return nothing

	# Advance column
	next = (i, j + 1)

	value = cartesian_from_tril(TriLayout(e), i, j)
	return value, next
end


"""
$(TYPEDSIGNATURES)

Iterate over linear and 2d indices of the lower triangle of an `n` by `n` matrix
in tandem.

This should be equivalent to the generator expression:

    (lin2car(layout, i) for i in 1:nelems(n))
"""
tri_indices(layout::TriLayout, n::Int) = TriIndexIterator{typeof(layout)}(n)
tri_indices(n::Int) = tri_indices(TriLower{true}(), n)
