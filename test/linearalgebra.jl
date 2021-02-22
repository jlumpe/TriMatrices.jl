using TriMatrices
using TriMatrices: wraptri
using TriMatrices.Indexing: trinum
using LinearAlgebra


const _N = 4
const _data = 1:trinum(_N)
const m_upper = TriMatrix(TriUpper(), _N, _data)
const m_lower = TriMatrix(TriLower(), _N, _data)
const m_sym = TriMatrix(TriSymmetric(), _N, _data)

const _data_nd = 1:trinum(_N - 1)
const _diag = -1
const m_upper_nd = TriMatrix(TriUpper{false}(), _N, _data_nd; diag=_diag)
const m_lower_nd = TriMatrix(TriLower{false}(), _N, _data_nd; diag=_diag)
const m_sym_nd = TriMatrix(TriSymmetric{false}(), _N, _data_nd; diag=_diag)


@testset "Transpose" begin
	t_upper = transpose(m_upper)
	@test t_upper isa TriMatrix{Int, TriLower{true}}
	@test parent(t_upper) === parent(m_upper)

	t_lower = transpose(m_lower)
	@test t_lower isa TriMatrix{Int, TriUpper{true}}
	@test parent(t_lower) === parent(m_lower)

	@test transpose(m_sym) === m_sym

	t_upper_nd = transpose(m_upper_nd)
	@test t_upper_nd isa TriMatrix{Int, TriLower{false}}
	@test parent(t_upper_nd) === parent(m_upper_nd)
	@test t_upper_nd.diag == m_upper_nd.diag

	t_lower_nd = transpose(m_lower_nd)
	@test t_lower_nd isa TriMatrix{Int, TriUpper{false}}
	@test parent(t_lower_nd) === parent(m_lower_nd)
	@test t_lower_nd.diag == m_lower_nd.diag

	@test transpose(m_sym_nd) === m_sym_nd
end


@testset "wraptri" begin
	cases = [
		[m_upper, m_upper_nd] => UpperTriangular,
		[m_lower, m_lower_nd] => LowerTriangular,
		[m_sym, m_sym_nd] => Symmetric,
	]

	for (ms, T) in cases
		for m in ms
			wrapped = wraptri(m)
			@test wrapped isa T
			@test parent(wrapped) === m
		end
	end
end