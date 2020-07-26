using TriMatrices
using TriMatrices: hasdiag
using TriMatrices.Testing

const N_VALS = [0, 1, 2, 5, 25]
const T_VALS = [Float64, Int]
# const DIAG_VALS = [1, 10]


@testset "Construct uninitialized" begin
	for L in LAYOUT_TYPES_P
		for n in N_VALS
			for T in T_VALS
				m = TriMatrix{L, T}(undef, n)
				@test typeof(m) === TriMatrix{L, T, Vector{T}}
				@test size(m) == (n, n)
				@test m.diag === zero(T)

				m = TriMatrix{L, T}(undef, n, 1)
				@test m.diag === one(T)
			end

			# Type defaults to Float64
			@test typeof(TriMatrix{L}(undef, n)) == TriMatrix{L, Float64, Vector{Float64}}
		end
	end
end


@testset "similar" begin
	for L in LAYOUT_TYPES_P
		for n in N_VALS
			for T in T_VALS
				m = TriMatrix{L, T}(undef, n, 1)
				m2 = similar(m)
				@test typeof(m2) == typeof(m)
				@test size(m2) == (n, n)
				@test m2.diag == m.diag

				# Change type
				T2 = Float32
				m3 = similar(m, T2)
				@test typeof(m3) == TriMatrix{L, T2, Vector{T2}}
				@test size(m3) == (n, n)
				@test m3.diag == m.diag

				# Change diag
				m4 = similar(m, diag=2)
				@test typeof(m4) == typeof(m)
				@test size(m4) == (n, n)
				@test m4.diag == 2
			end
		end
	end
end


@testset "Construct from AbstractMatrix" begin
	for L in LAYOUT_TYPES_P
		for n in N_VALS
			for T in T_VALS
				data, m = make_test_matrix(L(), n, T=T)

				m2 = TriMatrix{L}(m)
				@test m2 isa TriMatrix{L, T}
				@test m2 == m

				T2 = Float32
				m3 = TriMatrix{L}(m, T2)
				@test m3 isa TriMatrix{L, T2}
				@test m3 == m
			end
		end
	end
end


@testset "Construct filled" begin
	for L in LAYOUT_TYPES_P
		for n in N_VALS
			for T in T_VALS
				m = zeros(T, L(), n)
				@test m isa TriMatrix{L, T} && size(m) == (n, n)
				@test all(m.data .== 0)

				m = ones(T, L(), n)
				@test m isa TriMatrix{L, T} && size(m) == (n, n)
				@test all(m.data .== 1)

				m = fill(T(3), L(), n)
				@test m isa TriMatrix{L, T} && size(m) == (n, n)
				@test all(m.data .== 3)
			end
		end

		@test zeros(L(), 10) isa TriMatrix{L, Float64}
		@test ones(L(), 10) isa TriMatrix{L, Float64}
	end
end


@testset "Basic properties" begin
	# TODO
end


@testset "getindex and setindex!" begin
	# TODO
end


@testset "LinearAlgebra" begin
	# TODO
end
