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
				m = TriMatrix{T}(L(), undef, n)
				@test typeof(m) === TriMatrix{T, L, Vector{T}}
				@test size(m) == (n, n)
				@test m.diag === zero(T)

				m = TriMatrix{T}(L(), undef, n, 1)
				@test m.diag === one(T)
			end

			# Type defaults to Float64
			@test typeof(TriMatrix(L(), undef, n)) == TriMatrix{Float64, L, Vector{Float64}}
		end
	end
end


@testset "Construct from AbstractMatrix" begin
	for L in LAYOUT_TYPES_P
		for n in N_VALS
			for T in T_VALS
				data, m = make_test_matrix(L(), n, T=T)

				m2 = TriMatrix(L(), m)
				@test m2 isa TriMatrix{T, L}
				@test m2 == m

				T2 = Float32
				m3 = TriMatrix(L(), m, T2)
				@test m3 isa TriMatrix{T2, L}
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
				@test m isa TriMatrix{T, L} && size(m) == (n, n)
				@test all(m.data .== 0)

				m = ones(T, L(), n)
				@test m isa TriMatrix{T, L} && size(m) == (n, n)
				@test all(m.data .== 1)

				m = fill(T(3), L(), n)
				@test m isa TriMatrix{T, L} && size(m) == (n, n)
				@test all(m.data .== 3)
			end
		end

		@test zeros(L(), 10) isa TriMatrix{Float64, L}
		@test ones(L(), 10) isa TriMatrix{Float64, L}
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
