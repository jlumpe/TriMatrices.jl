using TriMatrices
using TriMatrices.Testing


@testset "make_test_matrix" begin
	@test make_test_matrix(TriLower{true}(), 3) == (
		[10, 20, 30, 40, 50, 60],
		[10  0  0;
		 20 30  0;
		 40 50 60]
	)
	@test make_test_matrix(TriLower{false}(), 3) == (
		[10, 20, 30],
		[-1  0  0;
		 10 -1  0;
		 20 30 -1]
	)

	@test make_test_matrix(TriUpper{true}(), 3) == (
		[10, 20, 30, 40, 50, 60],
		[10 20 40;
		 0 30 50;
		 0  0 60]
	)
	@test make_test_matrix(TriUpper{false}(), 3) == (
		[10, 20, 30],
		[-1 10 20;
		 0 -1 30;
		 0  0 -1]
	)

	@test make_test_matrix(TriSymmetric{true}(), 3) == (
		[10, 20, 30, 40, 50, 60],
		[10 20 40;
		 20 30 50;
		 40 50 60]
	)
	@test make_test_matrix(TriSymmetric{false}(), 3) == (
		[10, 20, 30],
		[-1 10 20;
		 10 -1 30;
		 20 30 -1]
	)
end
