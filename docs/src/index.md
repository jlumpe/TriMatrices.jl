# TriMatrices.jl


[`TriMatrix`](@ref) is a subtype of `AbstractMatrix` which stores only the values on one side of the
diagonal, thus reducing the memory required by about half.
The data is stored in an underlying 1D array.
It is intended for large matrices where the amount of memory saved is significant.


## Layouts

The `TriMatrix` type is parameterized by subtype of [`TriMatrices.TriLayout`](@ref) which determines how the
elements of the underlying data array are mapped to the 2D matrix. This is one of
[`TriUpper`](@ref), [`TriLower`](@ref), or [`TriSymmetric`](@ref).

This is illustrated by using a unit range for the data array:

```julia-repl
julia> n = 4

julia> data = 1:10

julia> TriMatrix(TriUpper(), n, data)
4×4 TriMatrix{Int64,TriUpper{true},UnitRange{Int64}}:
 1  2  4   7
 ⋅  3  5   8
 ⋅  ⋅  6   9
 ⋅  ⋅  ⋅  10

julia> TriMatrix(TriLower(), n, data)
4×4 TriMatrix{Int64,TriLower{true},UnitRange{Int64}}:
 1  ⋅  ⋅   ⋅
 2  3  ⋅   ⋅
 4  5  6   ⋅
 7  8  9  10

julia> TriMatrix(TriSymmetric(), n, data)
4×4 TriMatrix{Int64,TriSymmetric{true},UnitRange{Int64}}:
 1  2  4   7
 2  3  5   8
 4  5  6   9
 7  8  9  10
```


## Creating TriMatrix instances

In the below examples, the size of the resulting matrix is `n` by `n`.

* Wraping an existing data array (as shown above): `TriMatrix(layout::TriLayout, n, data::AbstractVector)`
* Uninitialized with eltype `T`: `TriMatrix{T=Float64}(layout::TriLayout, undef, n)`
* Copied from an existing square matrix `m`: `TriMatrix{T=eltype(m)}(layout::TriLayout, m::AbstractMatrix)`
* Filled with ones: `ones(T::Type=Float64, layout::TriLayout, n)`, same for `zeros`
* Filled with the value `x`: `fill(x, layout::TriLayout, n)`


## Matrices with constant diagonal

`TriLayout` subtypes have a single boolean type parameter which indicates whether elements along
the diagonal are explicitly stored along with the rest (`true`, the default) or filled with a
constant value (`false`). For the 2nd case, constructor methods take a `diag` keyword argument
to set this value:

```julia-repl
julia> TriMatrix(TriUpper{false}(), 4, 1:6; diag=-1)
4×4 TriMatrix{Int64,TriUpper{false},UnitRange{Int64}}:
 -1   1   2   4
  ⋅  -1   3   5
  ⋅   ⋅  -1   6
  ⋅   ⋅   ⋅  -1
```

This can be used with `TriSymmetric{false}` and `diag=0` or `diag=1` for a distance or
correlation matrix, respectively, or with `TriUpper{false}`/`TriLower{false}` and
`diag=1` for a unit triangular matrix.
The `diag` argument can also be passed to functions like `ones`, `zeros`, and `fill`.
It defaults to zero.
