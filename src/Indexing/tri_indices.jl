export tri_indices


"""
Iterator over Cartesian indices of the stored region of a `TriMatrix` in an
order corresponding to the elements of its data array. Created with the
[`tri_indices`](@ref) function.
"""
struct TriIndexIterator{L<:TriLayout}
	n::Int

	TriIndexIterator(::L, n::Int) where {L<:TriLayout} = new{L}(n)
end

Base.IteratorSize(::Type{<:TriIndexIterator}) = Base.HasLength()
Base.length(ti::TriIndexIterator{L}) where L = nelems(L, ti.n)
Base.IteratorEltype(::Type{<:TriIndexIterator}) = Base.HasEltype()
Base.eltype(::Type{<:TriIndexIterator}) = Tuple{Int, Int}
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

	value = cartesian_from_tril(L(), i, j)
	return value, next
end


"""
$(TYPEDSIGNATURES)

Iterate over the cartesian indices of a `TriMatrix` corresponding to the linear
indices of its data array, in order.

This should be equivalent to the following generator expression:

    (lin2car(layout, i) for i in 1:nelems(layout, n))
"""
tri_indices(layout::TriLayout, n::Integer) = TriIndexIterator(layout, n)
