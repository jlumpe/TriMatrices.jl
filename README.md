# TriMatrices.jl

[![Build Status](https://github.com/jlumpe/TriMatrices.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/jlumpe/TriMatrices.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jlumpe.github.io/TriMatrices.jl/stable)


Julia package for storing large triangular or symmetric matrices in non-redundant format.

`TriMatrix` is a subtype of `AbstractMatrix` which stores only the values on one side of the
diagonal, thus reducing the memory required by about half.
The data is stored in an underlying 1D array.
It is intended for large matrices where the amount of memory saved is significant.


## Example

```julia-repl
julia> A = rand(1000, 1000);

julia> A += A'  # Make symmetric;

julia> Base.summarysize(A)  # Number of bytes A occupies in memory - about 8MB
8000040

julia> using TriMatrices

julia> A2 = TriMatrix(TriSymmetric(), A);

julia> summary(stdout, A2)
100×100 TriMatrix{Int64,TriSymmetric{true},Array{Int64,1}}

julia> A == A2  # Arrays contain the same elements
true

juilia> Base.summarysize(A2)  # About 4MB
4004064
```
