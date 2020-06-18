# Useful integer math stuff

export trinum, triinv, triinv_rem


"""
$(TYPEDSIGNATURES)

Get the integer value of `floor(sqrt(n))`.
"""
@inline isqrt(n::Integer) = floor(typeof(n), sqrt(n))


"""
$(TYPEDSIGNATURES)

Get the nth triangular number.
"""
@inline trinum(n::Integer) = n * (n + 1) ÷ 2


"""
    triinv(t::Integer, strict=false)

Get `n` such that `t` is the `n`th triangular number.
If `strict=true` and `t` is not a triangular number, throw a `DomainError`.
"""
triinv(t::Integer, strict::Bool=false) = triinv(t, Val(strict))

function triinv(t::Integer, ::Val{strict}) where strict
	s = 8 * t + 1
	rs = isqrt(s)
	strict && rs^2 != s && throw(DomainError(t, "Not a triangular number"))
	return (rs - 1) ÷ 2
end


"""
$(TYPEDSIGNATURES)

Get `n` and `r >= 0` such that `t == trinum(n) + r`.
"""
function triinv_rem(t::Integer)
	n = triinv(t, Val{false}())
	return n, t - trinum(n)
end
