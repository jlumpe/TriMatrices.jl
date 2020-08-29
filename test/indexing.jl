using TriMatrices
using TriMatrices.Indexing
using TriMatrices.Testing


@testset "tri and triinv" begin
	for n in [0:100; 2 .^ 10:10:30]
		t = trinum(n)

		@test t == sum(1:n)

		@test triinv(t) == n
		@test triinv(t) == n
		@test triinv_strict(t) == n
	end

	for t in [0:100; 2 .^ 10:10:30]
		n, r = triinv_rem(t)
		@test n >= 0 && (0 <= r <= n)
		@test trinum(n) + r == t
		@test triinv(t) == n
		r == 0 || @test_throws DomainError triinv_strict(t)
	end
end


@testset "Index checking and conversion" begin
	for L in LAYOUT_TYPES_P
		layout = L()
		for n in 0:10
			data, tmat = make_test_matrix(layout, n)

			# By cartesian index
			for r in 1:n, c in 1:n
				v = tmat[r, c]

				if v > 0  # Element is stored in data array
					@test check_tri_index(Bool, layout, r, c)
					check_tri_index(layout, r, c)

					i = car2lin(layout, r, c)
					@test car2lin_unchecked(layout, r, c) == i

					@test tmat[r, c] == data[i]

				else  # Not stored
					@test !check_tri_index(Bool, layout, r, c)
					@test_throws DomainError check_tri_index(layout, r, c)
					@test_throws DomainError car2lin(layout, r, c)
				end
			end

			# By data index
			for i in eachindex(data)
				r, c = lin2car(layout, i)
				@test data[i] == tmat[r, c]
				@test car2lin(layout, r, c) == i
			end
		end
	end
end


@testset "tri_indices" begin
	for L in LAYOUT_TYPES_P
		layout = L()
		for n in 0:10
			data, tmat = make_test_matrix(layout, n)

			ti = tri_indices(layout, n)
			@test length(ti) == length(data)

			i = 0
			for (r, c) in tri_indices(layout, n)
				i += 1
				@test data[i] == tmat[r, c]
			end
			@test i == length(data)
		end
	end
end
