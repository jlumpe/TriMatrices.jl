# Useful integer math stuff

export trinum, triinv, triinv_strict, triinv_rem


"""
$(TYPEDSIGNATURES)

Get the integer value of `floor(sqrt(n))`.
"""
@inline isqrt(n::Integer) = floor(typeof(n), sqrt(n))


"""
$(TYPEDSIGNATURES)

Get the `n`th triangular number.

# See also
[`triinv`](@ref)
"""
@inline trinum(n::Integer) = n * (n + 1) ÷ 2


"""
$(TYPEDSIGNATURES)

Get the largest `n` such that `t >= trinum(n)`.

# See also
[`triinv_strict`](@ref), [`trinum`](@ref)
"""
@inline triinv(t::Integer) = _triinv(t, Val{false}())


"""
$(TYPEDSIGNATURES)

Get `n` such that `t` is the `n`th triangular number.
If `t` is not a triangular number, throw a `DomainError`.

# See also
[`triinv`](@ref), [`trinum`](@ref)
"""
@inline triinv_strict(t::Integer) = _triinv(t, Val{true}())

function _triinv(t::Integer, ::Val{strict}) where strict
	t >= 0 || throw(DomainError(t, "Must be nonnegative"))
	s = 8 * t + 1
	rs = isqrt(s)
	strict && rs^2 != s && throw(DomainError(t, "Not a triangular number"))
	return (rs - 1) ÷ 2
end


"""
$(TYPEDSIGNATURES)

Get `n` and `r >= 0` such that `t == trinum(n) + r` and `t < trinum(n + 1)`.

# See also
[`triinv`](@ref), [`trinum`](@ref)
"""
function triinv_rem(t::Integer)
	n = triinv(t)
	return n, t - trinum(n)
end
