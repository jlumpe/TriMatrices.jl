"""Helper code for running tests."""
module Testing

using LinearAlgebra: diagind
using TriMatrices
using TriMatrices: TriLayout, nelems, hasdiag


export LAYOUT_TYPES, LAYOUT_TYPES_P, make_test_matrix


const LAYOUT_TYPES = [TriLower, TriUpper, TriSymmetric]
const LAYOUT_TYPES_P = [L{D} for D in [true, false] for L in LAYOUT_TYPES]


"""
Make an `n` by `n` `Array` to test TriLayout and TriMatrix methods.

Values are as follows:
* Elements which are stored in matrix data have value `datafunc(i)` where `i`
  is the linear index in the data array. The mapping function is to make sure
  we're not confusing the value with its index.
* Elements along the diagonal have value `diag` if the diagonal is not stored
  for the given layout.
* All other elements are zero.

Returns `(data, matrix)`.
"""
function make_test_matrix(layout::TriLayout, n; diag=-1, datafunc=i -> 10i)
    data = map(datafunc, 1:nelems(layout, n))

    matrix = zeros(Int, n, n)
    matrix[diagind(matrix)] .= diag

    k = 0
    for i in 1:n, j in 1:i
        !hasdiag(layout) && i == j && continue
        k += 1
        layout isa TriUpper || (matrix[i, j] = data[k])
        layout isa TriLower || (matrix[j, i] = data[k])
    end
    @assert k == length(data)

    return data, matrix
end


end  # module
