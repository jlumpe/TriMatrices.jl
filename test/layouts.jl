using TriMatrices
using TriMatrices: trinum, hasdiag, nelems


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
