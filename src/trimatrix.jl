"""
$(TYPEDEF)

A triangular or symmetric matrix which stores data non-redundantly in a
contiguous linear array.
"""
struct TriMatrix{T, L<:TriLayout, A} <: AbstractMatrix{T}
	n::Int
	data::A
	diag::T

	"""
		TriMatrix(layout::TriLayout, n::Integer, data::AbstractVector{T}; diag=zero(T))

	Create an `n` by `n` `TriMatrix` using the existing 1D array `data` for storage.

	The length of `data` must match `n`, use [`TriMatrices.nelems`](@ref) to determine this value.
	"""
	function TriMatrix(layout::TriLayout, n::Integer, data::AbstractVector; diag=zero(eltype(data)))
		Base.require_one_based_indexing(data)
		length(data) == nelems(layout, n) || error("Matrix size does not match size of data array")
		return new{eltype(data), typeof(layout), typeof(data)}(n, data, diag)
	end
end


const TriUpperMatrix{T, D, A} = TriMatrix{T, TriUpper{D}, A}
const TriLowerMatrix{T, D, A} = TriMatrix{T, TriLower{D}, A}
const TriULMatrix{T, D, A} = TriMatrix{T, <:TriUL{D}, A}
const TriSymmetricMatrix{T, D, A} = TriMatrix{T, TriSymmetric{D}, A}



########################################
# Constructors
########################################


# Construct uninitialized
"""
	TriMatrix{T}(layout::TriLayout, undef, n::Integer; diag=zero(T))
	TriMatrix(layout::TriLayout, undef, n::Integer; diag=0.)

Create an uninitialized `n` by `n` `TriMatrix`.

`T` defaults to `Float64` if not specified.
"""
function TriMatrix{T}(layout::TriLayout, ::UndefInitializer, n::Integer; diag=zero(T)) where T
	data = Array{T}(undef, nelems(layout, n))
	return TriMatrix(layout, n, data, diag=convert(T, diag))
end

TriMatrix(layout::TriLayout, ::UndefInitializer, n::Integer; kw...) = TriMatrix{Float64}(layout, undef, n; kw...)


"""
	TriMatrix{T}(layout::TriLayout, m::AbstractMatrix; diag=zero(T))
	TriMatrix(layout::TriLayout, m::AbstractMatrix; [diag])

Create a `TriMatrix` by copying data from a square matrix `m`.

If unspecified `T` defaults to `eltype(m)`.
"""
function TriMatrix{T}(layout::TriLayout, m::AbstractMatrix;
                      diag=isempty(m) ? zero(T) : convert(T, m[1, 1])) where T
	n = LinearAlgebra.checksquare(m)

	tm = TriMatrix{T}(layout, undef, n, diag=diag)
	copyto!(tm, m)

	return tm
end

TriMatrix(layout::TriLayout, m::AbstractMatrix{T}; kw...) where T = TriMatrix{T}(layout, m; kw...)


"""
	fill(x, layout::TriLayout, n::Integer; diag=x)

Create an `n` by `n` [`TriMatrix`](@ref) with the stored data region filled with the value `x`.
"""
function Base.fill(x, layout::TriLayout, n::Integer; diag=x)
	mat = TriMatrix{typeof(x)}(layout, undef, n; diag=diag)
	fill!(mat.data, x)
	return mat
end


"""
	ones(T::Type=Float64, layout::TriLayout, n::Integer; diag=x)

Create an `n` by `n` [`TriMatrix`](@ref) with the stored data region filled with ones.
"""
Base.ones(T::Type, layout::TriLayout, n::Integer; diag=one(T)) = fill(one(T), layout, n, diag=diag)
Base.ones(layout::TriLayout, n::Integer; kw...) = ones(Float64, layout, n; kw...)


"""
	zeros(T::Type=Float64, layout::TriLayout, n::Integer; diag=x)

Create an `n` by `n` [`TriMatrix`](@ref) with the stored data region filled with zeros.
"""
Base.zeros(T::Type, layout::TriLayout, n::Integer; diag=zero(T)) = fill(zero(T), layout, n, diag=diag)
Base.zeros(layout::TriLayout, n::Integer; kw...) = zeros(Float64, layout, n; kw...)


# Construct similar
Base.similar(m::TriMatrix, new_eltype::Type, dims::Dims) = similar(m.data, new_eltype, dims)



########################################
# General methods
########################################

"""
	TriMatrices.TriLayout(m::TriMatrix)

Get the `TriLayout` of the given `TriMatrix` instance.
"""
TriLayout(m::TriMatrix) = TriLayout(typeof(m))
TriLayout(::Type{<:TriMatrix{T, L}}) where {T, L} = L()

TriMatrices.Indexing.tri_indices(m::TriMatrix) = tri_indices(TriLayout(m), m.n)


Base.size(m::TriMatrix) = (m.n, m.n)
Base.IndexStyle(::Type{<:TriMatrix}) = IndexCartesian()

"""
	parent(m::TriMatrix)

Get the underlying data array of the [`TriMatrix`](@ref) `m`.
"""
Base.parent(m::TriMatrix) = m.data


# Copy data from general matrix into TriMatrix
function Base.copyto!(dest::TriMatrix, src::AbstractMatrix)
	n = dest.n
	size(src) == (n, n) || error("Matrix sizes do not match")

	for (i, (r, c)) in enumerate(tri_indices(dest))
		@inbounds dest.data[i] = src[r, c]
	end

	return dest
end



########################################
# Indexing
########################################

"""
$(TYPEDSIGNATURES)

Get the value of an element within the stored data region of a
[`TriMatrix`](@ref). Should be faster than `Base.getindex` (even with
`@inbounds` applied) but results are undefined if `check_tri_index(mat, i, j)`
does not pass.
"""
function getindex_tri_unsafe(mat::TriMatrix, i::Integer, j::Integer)
	ii = car2lin_unchecked(TriLayout(mat), i, j)
	return @inbounds mat.data[ii]
end
getindex_tri_unsafe(mat::TriMatrix, idx::CartesianIndex{2}) = getindex_tri_unsafe(mat, idx[1], idx[2])

@Base.propagate_inbounds function Base.getindex(mat::TriMatrix{T,L},
                                                i::Integer, j::Integer,
                                                ) where {T, L}

	@boundscheck checkbounds(mat, i, j)

	!hasdiag(L) && i == j && return mat.diag  # Diagonal fill value
	L <: TriLower && i < j && return zero(T)  # Lower triangular above the diagonal
	L <: TriUpper && i > j && return zero(T)  # Upper triangular below the diagonal

	return getindex_tri_unsafe(mat, i, j)
end


"""
$(TYPEDSIGNATURES)

Set the value of an element within the stored data region of a
[`TriMatrix`](@ref). Should be faster than `Base.setindex!` (even with
`@inbounds` applied) but results are undefined if `check_tri_index(mat, i, j)`
does not pass.
"""
function setindex_tri_unsafe!(mat::TriMatrix, v, i::Integer, j::Integer)
	ii = car2lin_unchecked(TriLayout(mat), i, j)
	return @inbounds mat.data[ii] = v
end
setindex_tri_unsafe!(mat::TriMatrix, v, idx::CartesianIndex{2}) = setindex_tri_unsafe!(mat, v, idx[1], idx[2])

@Base.propagate_inbounds function Base.setindex!(mat::TriMatrix{T,L}, v,
                                                 i::Integer, j::Integer,
                                                ) where {T, L}

	@boundscheck checkbounds(mat, i, j)

	if !hasdiag(L) && i == j
		v != mat.diag && error("Cannot change value on diagonal for TriMatrix with layout $L")
		return mat.diag
	end

	if L <: TriLower && i < j
		v != 0 && error("Cannot assign non-zero value to index above diagonal for TriMatrix with layout TriLower")
		return zero(T)
	end

	if L <: TriUpper && i > j
		v != 0 && error("Cannot assign non-zero value to index below diagonal for TriMatrix with layout TriUpper")
		return zero(T)
	end

	return setindex_tri_unsafe!(mat, v, i, j)
end



########################################
# Display
########################################

# Renders entries on zero side of diagonal as dots, like triangular types in LinearAlgebra
function Base.replace_in_print_matrix(::TriMatrix{T, <:TriLower}, i::Integer, j::Integer, s::AbstractString) where T
	i >= j ? s : Base.replace_with_centered_mark(s)
end
function Base.replace_in_print_matrix(::TriMatrix{T, <:TriUpper}, i::Integer, j::Integer, s::AbstractString) where T
	i <= j ? s : Base.replace_with_centered_mark(s)
end
