# Convert between 1D and 2D indexing for triangular arrays

"""
$(TYPEDSIGNATURES)

Convert a matrix-style row/column index to linear index of lower triangular values.

Linear indices are arranged as follows:

    1
    2 3
    4 5 6
    ...

The domain of this function is `i >= 1`, `1 <= j <= i`.

This function is the inverse of [`tri2mat`](@ref).
"""
function mat2tri(i, j)
    (i >= 1 && (1 <= j <= i)) || throw(DomainError("Invalid lower triangular index: ($i, $j)"))
    return mat2tri_unchecked(i, j)
end

@inline mat2tri_unchecked(i, j) = tri(i - 1) + j


"""
$(TYPEDSIGNATURES)

Convert a linear index of lower triangular values to a row/column matrix index.

Linear indices are arranged as follows:

    1
    2 3
    4 5 6
    ...

This function is the inverse of [`mat2tri`](@ref).
"""
function tri2mat(n)
    n, r = triinv_rem(n - 1)
    return (n + 1, r + 1)
end
