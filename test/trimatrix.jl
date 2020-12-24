using TriMatrices
using TriMatrices: hasdiag, TriLayout, getindex_tri_unsafe, setindex_tri_unsafe!
using TriMatrices.Testing

const N_VALS = [0, 1, 2, 5, 25]
const T_VALS = [Float64, Int]
# const DIAG_VALS = [1, 10]


@testset "Construct uninitialized" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			for T in T_VALS
				m = TriMatrix{T}(layout, undef, n)
				@test typeof(m) === TriMatrix{T, L, Vector{T}}
				@test size(m) == (n, n)
				@test m.diag === zero(T)

				m = TriMatrix{T}(layout, undef, n, diag=1)
				@test m.diag === one(T)
			end

			# Type defaults to Float64
			@test typeof(TriMatrix(layout, undef, n)) == TriMatrix{Float64, L, Vector{Float64}}
		end
	end
end


@testset "Construct from AbstractMatrix" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			for T in T_VALS
				data, m = make_test_matrix(layout, n, T=T)

				m2 = TriMatrix(layout, m)
				@test m2 isa TriMatrix{T, L}
				@test m2 == m

				T2 = Float32
				m3 = TriMatrix(layout, m, T2)
				@test m3 isa TriMatrix{T2, L}
				@test m3 == m
			end
		end

		@test_throws DimensionMismatch TriMatrix(layout, zeros(3, 4))
	end
end


@testset "Construct filled" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			for T in T_VALS
				cases = [
					(zeros(T, layout, n), zero(T)),
					(ones(T, layout, n), one(T)),
					(fill(T(3), layout, n), T(3)),
				]

				for (m, v) in cases
					@test m isa TriMatrix{T, L}
					@test size(m) == (n, n)
					@test all(x -> x === v, m.data)
				end
			end
		end

		@test zeros(layout, 10) isa TriMatrix{Float64, L}
		@test ones(layout, 10) isa TriMatrix{Float64, L}
	end
end


@testset "Basic properties" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			for T in T_VALS
				m = TriMatrix{T}(layout, undef, n)

				@test size(m) == (n, n)
				@test eltype(m) == eltype(typeof(m)) == T
				@test Base.IndexStyle(m) == Base.IndexStyle(typeof(m)) == Base.IndexCartesian()
				@test Base.parent(m) === m.data

				@test TriLayout(m) == TriLayout(typeof(m)) == layout
			end
		end
	end
end


@testset "similar" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			for T in T_VALS
				am, tm = make_test_matrix_pair(layout, n, T=T)

				# Explicit calls to similar()
				sm = similar(tm)
				@test sm isa Array{T} && size(sm) == size(tm)

				sm = similar(tm, Float32)
				@test sm isa Array{Float32} && size(sm) == size(tm)

				sm = similar(tm, (10,))
				@test sm isa Array{T} && size(sm) == (10,)

				sm = similar(tm, Float32, (10,))
				@test sm isa Array{Float32} && size(sm) == (10,)

				# Implicit calls when allocating for result
				if n > 0
					sm = tm[1, :]
					@test sm isa Vector{T} && sm == am[1, :]
				end
			end
		end
	end
end


@testset "getindex" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			n == 0 && continue
			for T in T_VALS
				am, tm = make_test_matrix_pair(layout, n, T=T)

				for i in 1:n, j in 1:n
					@test tm[i, j] === am[i, j]
					@inferred tm[i, j]

					idx = CartesianIndex(i, j)
					@test tm[idx] === am[i, j]
					@inferred tm[idx]

					if check_tri_index(Bool, layout, i, j)
						@test getindex_tri_unsafe(tm, i, j) === am[i, j]
						@test getindex_tri_unsafe(tm, idx) === am[i, j]
					end
				end
			end
		end
	end
end


@testset "setindex!" begin
	for L in LAYOUT_TYPES_P
		layout = L()

		for n in N_VALS
			n == 0 && continue
			for T in T_VALS
				m = ones(layout, n, diag=-1)

				v = T(2)

				for i in 1:n, j in 1:n
					idx = CartesianIndex(i, j)

					if check_tri_index(Bool, layout, i, j)
						v += 1
						@test (m[i, j] = v) === v
						@test m[i, j] == v

						v += 1
						@test (m[idx] = v) === v
						@test m[i, j] == v

						v += 1
						@test setindex_tri_unsafe!(m, v, i, j) === v
						@test m[i, j] == v

						v += 1
						@test setindex_tri_unsafe!(m, v, idx) === v
						@test m[i, j] == v

					else
						# Setting to current value is OK
						@test (m[i, j] = m[i, j]) === m[i, j]
						# Changing is not
						@test_throws ErrorException m[i, j] += 1
					end
				end
			end
		end
	end
end


@testset "LinearAlgebra" begin
	# TODO
end
