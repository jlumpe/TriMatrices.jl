"""
$(TYPEDEF)

A triangular or symmetric matrix which stores data non-redundantly in contiguous
memory.
"""
struct TriMatrix{L<:TriLayout, T, A} <: AbstractMatrix{T}
	n::Int
	data::A
	diag::T

	function TriMatrix{L}(n::Int, data::A, diag=zero(T)) where {L<:TriLayout, T, A<:AbstractVector{T}}
		L isa DataType || throw(ArgumentError("Layout type parameter requires parameter"))
		Base.require_one_based_indexing(data)
		length(data) == nelems(L, n) || error("Matrix size does not match size of data array")
		return new{L, T, A}(n, data, diag)
	end
end


########################################
# Constructors
########################################

# May give layout as first argument instead of type parameter to any constructor
TriMatrix(::L, args...) where {L <: TriLayout} = TriMatrix{L}(args...)


# Construct uninitialized
function TriMatrix{L, T}(::UndefInitializer, n::Int, diag=zero(T)) where {L, T}
	data = Array{T}(undef, nelems(L, n))
	return TriMatrix{L}(n, data, diag)
end

function TriMatrix{L}(::UndefInitializer, n::Int, diag=0.) where {L}
	return TriMatrix{L, Float64}(undef, n, diag)
end

end


# Construct from existing matrix
function TriMatrix{L}(m::AbstractMatrix, T::Type=eltype(m), diag=isempty(m) ? zero(T) : convert(T, m[1, 1])) where L
	n = size(m, 1)
	size(m, 2) == n || error("Matrix is not square")

	tm = TriMatrix{L, T}(undef, n, diag)
	copyto!(tm, m)

	return tm
end


# Construct filled
function Base.fill(x, ::L, n::Int, diag=x) where {L <: TriLayout}
	mat = TriMatrix{L, typeof(x)}(undef, n, diag)
	fill!(mat.data, x)
	return mat
end

Base.ones(T::Type, layout::L, n, diag=one(T)) where {L <: TriLayout} = fill(one(T), layout, n, diag)
Base.ones(layout::L, n, diag=1.) where {L <: TriLayout} = ones(Float64, layout, n, diag)
Base.zeros(T::Type, layout::L, n, diag=zero(T)) where {L <: TriLayout} = fill(zero(T), layout, n, diag)
Base.zeros(layout::L, n, diag=0.) where {L <: TriLayout} = zeros(Float64, layout, n, diag)


# Construct similar
Base.similar(m::TriMatrix, new_eltype::Type, dims::Tuple) = similar(m.data, new_eltype, dims)
Base.similar(m::TriMatrix, dims::Tuple) = similar(m.data, dims)


########################################
# General methods
########################################

TriLayout(::Type{M}) where {L, M<:TriMatrix{L}} = L()
TriLayout(m::TriMatrix) = TriLayout(typeof(m))
TriMatrices.Indexing.tri_indices(m::TriMatrix) = tri_indices(TriLayout(m), m.n)


Base.size(m::TriMatrix) = (m.n, m.n)
Base.IndexStyle(::Type{<:TriMatrix}) = IndexCartesian()
Base.parent(m::TriMatrix) = m.data


@Base.propagate_inbounds function Base.getindex(mat::TriMatrix{L,T},
                                                r::Integer, c::Integer,
                                                ) where {L, T}

	@boundscheck checkbounds(mat, r, c)

	!hasdiag(L) && r == c && return mat.diag  # Diagonal fill value
	L <: TriLower && r < c && return zero(T)  # Lower triangular above the diagonal
	L <: TriUpper && r > c && return zero(T)  # Upper triangular below the diagonal

	i = car2lin_unchecked(L(), r, c)
	return mat.data[i]
end


@Base.propagate_inbounds function Base.setindex!(mat::TriMatrix{L,T},
                                                 v, r::Integer, c::Integer,
                                                 ) where {L, T}

	@boundscheck checkbounds(mat, r, c)

	!hasdiag(L) && r == c && v != mat.diag && error("Cannot change value on diagonal for TriMatrix with layout $L")

	if L <: TriLower && r < c
		v != 0 && error("Cannot assign non-zero value to index above diagonal for TriMatrix with layout TriLower")
		return zero(T)
	end

	if L <: TriUpper && c < r
		v != 0 && error("Cannot assign non-zero value to index below diagonal for TriMatrix with layout TriUpper")
		return zero(T)
	end

	i = car2lin_unchecked(L(), r, c)
	return mat.data[i] = v
end


# Copy data from general matrix into TriMatrix
function Base.copyto!(dest::TriMatrix{L}, src::AbstractMatrix) where {L}
	n = dest.n
	size(src) == (n, n) || error("Matrix sizes do not match")

	for (i, (r, c)) in enumerate(tri_indices(dest))
		dest.data[i] = src[r, c]
	end

	return dest
end


########################################
# Display
########################################

# Renders entries on zero side of diagonal as dots, like triangular types in LinearAlgebra
function Base.replace_in_print_matrix(::TriMatrix{<:TriLower}, i::Integer, j::Integer, s::AbstractString)
	i >= j ? s : Base.replace_with_centered_mark(s)
end
function Base.replace_in_print_matrix(::TriMatrix{<:TriUpper}, i::Integer, j::Integer, s::AbstractString)
	i <= j ? s : Base.replace_with_centered_mark(s)
end


########################################
# LinearAlgebra methods
########################################

"""
	$(FUNCTIONNAME)(m::TriMatrix)

Wrap a [`TriMatrix`](@ref) in the appropriate `LinearAlgebra` view type for
optimized performance.
"""
function wraptri end

wraptri(m::TriMatrix{<:TriLower}) = LowerTriangular(m)
wraptri(m::TriMatrix{<:TriUpper}) = UpperTriangular(m)
wraptri(m::TriMatrix{<:TriSymmetric}) = Symmetric(m)


Base.transpose(m::TriMatrix{TriUpper, L}) where L = TriMatrix{TriLower, L}(m.n, m.data, m.diag)
Base.transpose(m::TriMatrix{TriLower, L}) where L = TriMatrix{TriUpper, L}(m.n, m.data, m.diag)
Base.transpose(m::TriMatrix{TriSymmetric}) = m


LinearAlgebra.istriu(::TriMatrix{<:TriUpper}, n=0) = n >= 0 ? true : error("Not implemented")
LinearAlgebra.istril(::TriMatrix{<:TriLower}, n=0) = n <= 0 ? true : error("Not implemented")
LinearAlgebra.issymmetric(::TriMatrix{<:TriSymmetric}) = true


# TODO - override common matrix methods to operate on wraptri() value
# unary -
# binary +, -, *, /
# inv, dot, det, tr, eigen?
# mul!, lmul!, rmul!, ldiv!, rdiv!
