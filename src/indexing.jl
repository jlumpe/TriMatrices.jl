# Convert between 1D and 2D indexing for triangular arrays

"""
    check_tri_index(Bool, layout::TriLayout, i, j)::Bool

Check whether a row/column index is within the stored region of a
[`TriMatrix`](@ref) with the given layout.
"""
function check_tri_index end

@inline check_tri_index(::Type{Bool}, ::TriLower{true}, i, j) = j <= i
@inline check_tri_index(::Type{Bool}, ::TriLower{false}, i, j) = j < i
@inline check_tri_index(::Type{Bool}, ::TriUpper{true}, i, j) = i <= j
@inline check_tri_index(::Type{Bool}, ::TriUpper{false}, i, j) = i < j
@inline check_tri_index(::Type{Bool}, ::TriSymmetric{true}, i, j) = true
@inline check_tri_index(::Type{Bool}, ::TriSymmetric{false}, i, j) = i != j

"""
    check_tri_index(layout::TriLayout, i, j)

Check that a row/column index is within the stored region of a
[`TriMatrix`](@ref) with the given layout or throw a `DomainError`. This assumes
the index is otherwise valid for the corresponding matrix size.
"""
function check_tri_index(layout::TriLayout, i, j)
    check_tri_index(Bool, layout, i, j) || throw(DomainError("Invalid index for layout $(layout): ($i, $j)"))
end

# Convert cartesian index of given layout to cartesian index of TriLower{true}
# which has the same index in the data array.
# Assumes index is within the stored data region.
@inline cartesian_to_tril(::TriLower{true}, i, j) = (i, j)
@inline cartesian_to_tril(::TriLower{false}, i, j) = (i - 1, j)
@inline cartesian_to_tril(::TriUpper{true}, i, j) = (j, i)
@inline cartesian_to_tril(::TriUpper{false}, i, j) = (j - 1, i)
@inline function cartesian_to_tril(::TriSymmetric{D}, i, j) where D
    if i >= j
        return cartesian_to_tril(TriLower{D}(), i, j)
    else
        return cartesian_to_tril(TriUpper{D}(), i, j)
    end
end

# Inverse of cartesian_to_tril()
@inline cartesian_from_tril(::TriLower{true}, i, j) = (i, j)
@inline cartesian_from_tril(::TriLower{false}, i, j) = (i + 1, j)
@inline cartesian_from_tril(::TriUpper{true}, i, j) = (j, i)
@inline cartesian_from_tril(::TriUpper{false}, i, j) = (j, i + 1)
@inline cartesian_from_tril(::TriSymmetric{D}, i, j) where D = cartesian_from_tril(TriLower{D}(), i, j)


@inline car2lin_unchecked(i, j) = trinum(i - 1) + j
@inline car2lin_unchecked(layout::TriLayout, i, j) = car2lin_unchecked(cartesian_to_tril(layout, i, j)...)
# @inline car2lin_unchecked(::TriLower{true}, i, j) = car2lin_unchecked(i, j)
# @inline car2lin_unchecked(::TriLower{false}, i, j) = car2lin_unchecked(i - 1, j)
# @inline car2lin_unchecked(::TriUpper{D}, i, j) where D = car2lin_unchecked(TriLower{D}, j, i)
# function @inline car2lin_unchecked(::TriSymmetric{D}, i, j) where D
#     if i >= j
#         return car2lin_unchecked(TriLower{D}(), i, j)
#     else
#         return car2lin_unchecked(TriUpper{D}(), j, i)
#     end
# end

"""
$(TYPEDSIGNATURES)

Convert a Cartesian row/column index to a linear index of the data array of a
[`TriMatrix`](@ref) with the given layout.

This function is the inverse of [`lin2car`](@ref).
"""
function car2lin(layout::TriLayout, i, j)
    check_tri_index(layout, i, j)
    return car2lin_unchecked(layout, i, j)
end


"""
    $(FUNCTIONNAME)(layout::TriLayout, i)

Convert a linear index of the data array of a [`TriMatrix`](@ref) with the given
layout to the corresponding Cartesian index.

This function is the inverse of [`car2lin`](@ref).
"""
function lin2car(layout::TriLayout, i::Integer)
    r, c = lin2car(i)
    return cartesian_from_tril(layout, r, c)
end

function lin2car(i::Integer)
    n, i = triinv_rem(i - 1)
    return (n + 1, i + 1)
end


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
