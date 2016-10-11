



@inline split_exponent(q::Int64) = _split_exponent(q, UInt64(9), UInt64(31), UInt64(2))
@inline split_exponent(q::Int32) = _split_exponent(q, UInt32(6), UInt32(31), UInt32(2))
@inline split_exponent(q::Int16) = _split_exponent(q, UInt16(3), UInt16(31), UInt16(2))
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    m, q
end

"""
    ldexpk(a::FloatTypes, n::Int) -> Float

Computes `a × 2^n`.
"""
@inline function _ldexp{T<:FloatTypes}(x::T, k::Integer)
    q = k % inttype(T)
    bias = asint(T, exponent_bias(T))
    emax = asint(T, exponent_max(T))
    m, q = split_exponent(q)
    m += bias
    m = ifelse(m < 0, asint(T,0), m)
    m = ifelse(m > emax, emax, m)
    q += bias
    u = intasfloat(m)
    x = x*u*u*u*u
    u = intasfloat(q)
    return x*u
end
