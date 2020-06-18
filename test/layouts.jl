using TriMatrices
using TriMatrices: trinum, hasdiag, nelems


LAYOUT_TYPES = [TriLower, TriUpper, TriSymmetric]


for L in LAYOUT_TYPES
	n = 10

	# Shouldn't throw any errors
	L{true}()
	L{false}()

	# Default to true
	@test typeof(L()) === L{true}

	# Non-bool parameter
	@test_throws ArgumentError L{Int}()

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
