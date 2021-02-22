"""
$(TYPEDEF)

A triangular or symmetric matrix which stores data non-redundantly in a
contiguous linear array.
"""
struct TriMatrix{T, L<:TriLayout, A} <: AbstractMatrix{T}
	n::Int
	data::A
	diag::T

	function TriMatrix(::L, n::Int, data::A; diag=zero(T)) where {L<:TriLayout, T, A<:AbstractVector{T}}
		Base.require_one_based_indexing(data)
		length(data) == nelems(L, n) || error("Matrix size does not match size of data array")
		return new{T, L, A}(n, data, diag)
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
function TriMatrix{T}(layout::TriLayout, ::UndefInitializer, n::Int; diag=zero(T)) where T
	data = Array{T}(undef, nelems(layout, n))
	return TriMatrix(layout, n, data, diag=convert(T, diag))
end

TriMatrix(layout::TriLayout, ::UndefInitializer, n::Int; kw...) = TriMatrix{Float64}(layout, undef, n; kw...)


# Construct from existing matrix
function TriMatrix{T}(layout::TriLayout, m::AbstractMatrix;
                      diag=isempty(m) ? zero(T) : convert(T, m[1, 1])) where T
	n = LinearAlgebra.checksquare(m)

	tm = TriMatrix{T}(layout, undef, n, diag=diag)
	copyto!(tm, m)

	return tm
end

TriMatrix(layout::TriLayout, m::AbstractMatrix{T}; kw...) where T = TriMatrix{T}(layout, m; kw...)


# Construct filled
function Base.fill(x, layout::TriLayout, n::Int; diag=x)
	mat = TriMatrix{typeof(x)}(layout, undef, n; diag=diag)
	fill!(mat.data, x)
	return mat
end

Base.ones(T::Type, layout::TriLayout, n::Int; diag=one(T)) = fill(one(T), layout, n, diag=diag)
Base.ones(layout::TriLayout, n::Int; kw...) = ones(Float64, layout, n; kw...)
Base.zeros(T::Type, layout::TriLayout, n::Int; diag=zero(T)) = fill(zero(T), layout, n, diag=diag)
Base.zeros(layout::TriLayout, n::Int; kw...) = zeros(Float64, layout, n; kw...)


# Construct similar
Base.similar(m::TriMatrix, new_eltype::Type, dims::Dims) = similar(m.data, new_eltype, dims)



########################################
# General methods
########################################

TriLayout(::Type{<:TriMatrix{T, L}}) where {T, L} = L()
TriLayout(m::TriMatrix) = TriLayout(typeof(m))
TriMatrices.Indexing.tri_indices(m::TriMatrix) = tri_indices(TriLayout(m), m.n)


Base.size(m::TriMatrix) = (m.n, m.n)
Base.IndexStyle(::Type{<:TriMatrix}) = IndexCartesian()
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
