# API


## TriLayout

```@docs
TriMatrices.TriLayout
TriUpper
TriLower
TriSymmetric
TriMatrices.hasdiag
TriMatrices.nelems
TriMatrices.transpose_layout
```


## TriMatrix

```@docs
TriMatrix
TriMatrix(::TriMatrices.TriLayout, ::Integer, ::AbstractVector)
TriMatrix(::TriMatrices.TriLayout, ::UndefInitializer, ::Integer)
TriMatrix(::TriMatrices.TriLayout, ::AbstractMatrix)
Base.fill(::Any, ::TriMatrices.TriLayout, ::Integer)
Base.ones(::Type, ::TriMatrices.TriLayout, ::Integer)
Base.zeros(::Type, ::TriMatrices.TriLayout, ::Integer)
TriMatrices.TriLayout(::TriMatrix)
Base.parent(::TriMatrix)
TriMatrices.getindex_tri_unsafe
TriMatrices.setindex_tri_unsafe!
```


## Indexing

```@docs
TriMatrices.Indexing.check_tri_index
TriMatrices.Indexing.car2lin
TriMatrices.Indexing.car2lin_unchecked
TriMatrices.Indexing.lin2car
TriMatrices.Indexing.tri_indices
TriMatrices.Indexing.TriIndexIterator
TriMatrices.Indexing.trinum
TriMatrices.Indexing.triinv
TriMatrices.Indexing.triinv_strict
TriMatrices.Indexing.triinv_rem
```
