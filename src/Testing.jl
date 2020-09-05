"""Helper code for running tests."""
module Testing

using LinearAlgebra: diagind
using TriMatrices
using TriMatrices: TriLayout, nelems, hasdiag
using TriMatrices.Indexing: check_tri_index, car2lin


export LAYOUT_TYPES, LAYOUT_TYPES_P, make_test_matrix, make_test_matrix_pair


const LAYOUT_TYPES = [TriLower, TriUpper, TriSymmetric]
const LAYOUT_TYPES_P = [L{D} for D in [true, false] for L in LAYOUT_TYPES]


default_datafunc(i) = 10i

"""
Make an `n` by `n` `Array` to test TriLayout and TriMatrix methods.

First creates a 1D data vector of the appropriate size for the given layout and
value of `n` and files with values `datafunc(i)` for `i` in `1:nelems(layout, n)`.
The mapping function is to make sure we're not confusing the value with its index.

Then creates an `n` by `n` matrix with the folling values:
* Elements which are in the stored data region for the given `layout` have the
  value of the corresponding element of the data array.
* Elements along the diagonal have value `diag` if the diagonal is not stored
  for the given layout.
* All other elements are zero.

Returns `(data, matrix)`, where `matrix` is the `n` by `n` `Array` and `data`
is the data array for the corresponding `TriMatrix`.
"""
function make_test_matrix(layout::TriLayout, n::Int; diag=-1, datafunc=default_datafunc, T::Type=Int)
	data = collect(T, datafunc(i) for i in 1:nelems(layout, n))

	matrix = zeros(T, n, n)

	for i in 1:n, j in 1:n
		if check_tri_index(Bool, layout, i, j)
			matrix[i, j] = data[car2lin(layout, i, j)]
		elseif i == j
			matrix[i, j] = diag
		end
	end

	return data, matrix
end


"""
Create an `Array` using `make_test_matrix` along with a `TriMatrix` that should
contain the same values.

Returns a `(::Matrix, ::TriMatrix)` tuple.
"""
function make_test_matrix_pair(layout::TriLayout, n::Int; diag=-1, kw...)
	data, arraymat = make_test_matrix(layout, n; diag=diag, kw...)
	trimat = TriMatrix(layout, n, data, diag=diag)
	return arraymat, trimat
end


end  # module
