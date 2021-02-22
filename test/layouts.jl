using TriMatrices
using TriMatrices: trinum, hasdiag, nelems, transpose_layout


LAYOUT_TYPES = [TriLower, TriUpper, TriSymmetric]


@testset "Type parameter" begin
	for L in LAYOUT_TYPES
		# Shouldn't throw any errors
		L{true}()
		L{false}()

		# Default to true
		@test typeof(L()) === L{true}

		# Non-bool parameter
		@test_throws ArgumentError L{Int}()
	end
end


@testset "Methods and type traits" begin
	n = 10

	for L in LAYOUT_TYPES
		for D in [true, false]
			layout = L{D}()

			@test hasdiag(layout) === D
			@test hasdiag(typeof(layout)) === D

			ne = trinum(D ? n : n - 1)
			@test nelems(layout, n) == ne
			@test nelems(typeof(layout), n) == ne
		end

		# Type parameter needed
		@test_throws MethodError hasdiag(L)
		@test_throws MethodError nelems(L, n)
	end
end


@testset "transpose_layout" begin
	@test transpose_layout(TriLower) === TriUpper
	@test transpose_layout(TriUpper) === TriLower
	@test transpose_layout(TriSymmetric) === TriSymmetric
	@test transpose_layout(TriLower{true}) === TriUpper{true}
	@test transpose_layout(TriUpper{true}) === TriLower{true}
	@test transpose_layout(TriSymmetric{true}) === TriSymmetric{true}
	@test transpose_layout(TriLower{false}) === TriUpper{false}
	@test transpose_layout(TriUpper{false}) === TriLower{false}
	@test transpose_layout(TriSymmetric{false}) === TriSymmetric{false}

	# Check transpose of instance matches transpose of its type
	for L in LAYOUT_TYPES
		for D in [true, false]
			layout = L{D}()
			@test typeof(transpose_layout(layout)) === transpose_layout(typeof(layout))
		end
	end
end
