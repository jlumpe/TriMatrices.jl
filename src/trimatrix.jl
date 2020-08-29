"""
$(TYPEDEF)

A triangular or symmetric matrix which stores data non-redundantly in a
contiguous linear array.
"""
struct TriMatrix{T, L<:TriLayout, A} <: AbstractMatrix{T}
	n::Int
	data::A
	diag::T

	function TriMatrix(::L, n::Int, data::A, diag=zero(T)) where {L<:TriLayout, T, A<:AbstractVector{T}}
		Base.require_one_based_indexing(data)
		length(data) == nelems(L, n) || error("Matrix size does not match size of data array")
		return new{T, L, A}(n, data, diag)
	end
end


########################################
# Constructors
########################################


# Construct uninitialized
function TriMatrix{T}(layout::TriLayout, ::UndefInitializer, n::Int, diag=zero(T)) where T
	data = Array{T}(undef, nelems(layout, n))
	return TriMatrix(layout, n, data, convert(T, diag))
end

function TriMatrix(layout::TriLayout, ::UndefInitializer, n::Int, diag=0.)
	return TriMatrix{Float64}(layout, undef, n, diag)
end


# Construct from existing matrix
function TriMatrix(layout::TriLayout, m::AbstractMatrix, T::Type=eltype(m);
                   diag=isempty(m) ? zero(T) : convert(T, m[1, 1]))
	n = size(m, 1)
	size(m, 2) == n || error("Matrix is not square")

	tm = TriMatrix{T}(layout, undef, n, diag)
	copyto!(tm, m)

	return tm
end


# Construct filled
function Base.fill(x, layout::TriLayout, n::Int, diag=x)
	mat = TriMatrix{typeof(x)}(layout, undef, n, diag)
	fill!(mat.data, x)
	return mat
end

Base.ones(T::Type, layout::TriLayout, n, diag=one(T)) = fill(one(T), layout, n, diag)
Base.ones(layout::TriLayout, n, diag=1.) = ones(Float64, layout, n, diag)
Base.zeros(T::Type, layout::TriLayout, n, diag=zero(T)) = fill(zero(T), layout, n, diag)
Base.zeros(layout::TriLayout, n, diag=0.) = zeros(Float64, layout, n, diag)


# Construct similar
Base.similar(m::TriMatrix, new_eltype::Type, dims::Tuple) = similar(m.data, new_eltype, dims)
Base.similar(m::TriMatrix, dims::Tuple) = similar(m.data, dims)


########################################
# General methods
########################################

TriLayout(::Type{<:TriMatrix{T, L}}) where {T, L} = L()
TriLayout(m::TriMatrix) = TriLayout(typeof(m))
TriMatrices.Indexing.tri_indices(m::TriMatrix) = tri_indices(TriLayout(m), m.n)


Base.size(m::TriMatrix) = (m.n, m.n)
Base.IndexStyle(::Type{<:TriMatrix}) = IndexCartesian()
Base.parent(m::TriMatrix) = m.data


@Base.propagate_inbounds function Base.getindex(mat::TriMatrix{T,L},
                                                r::Integer, c::Integer,
                                                ) where {T, L}

	@boundscheck checkbounds(mat, r, c)

	!hasdiag(L) && r == c && return mat.diag  # Diagonal fill value
	L <: TriLower && r < c && return zero(T)  # Lower triangular above the diagonal
	L <: TriUpper && r > c && return zero(T)  # Upper triangular below the diagonal

	i = car2lin_unchecked(TriLayout(mat), r, c)
	return mat.data[i]
end


@Base.propagate_inbounds function Base.setindex!(mat::TriMatrix{T,L},
                                                 v, r::Integer, c::Integer,
                                                 ) where {T, L}

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

	i = car2lin_unchecked(TriLayout(mat), r, c)
	return mat.data[i] = v
end


# Copy data from general matrix into TriMatrix
function Base.copyto!(dest::TriMatrix, src::AbstractMatrix)
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
function Base.replace_in_print_matrix(::TriMatrix{T, <:TriLower}, i::Integer, j::Integer, s::AbstractString) where T
	i >= j ? s : Base.replace_with_centered_mark(s)
end
function Base.replace_in_print_matrix(::TriMatrix{T, <:TriUpper}, i::Integer, j::Integer, s::AbstractString) where T
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
