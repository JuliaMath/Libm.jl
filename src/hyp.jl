# overflow argument value
_over_sinhcosh(::Type{Float64}) = 710.0
_over_sinhcosh(::Type{Float32}) = 89f0

function sinh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    d = ddsub(d, ddrec(d))
    y = T(d) * T(0.5)
    y = abs(x) > _over_sinhcosh(T) ? typemax(T) : y
    y = isnan(y) ? typemax(T) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN)     : y
    return y
end

function cosh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    d = ddadd(d, ddrec(d))
    y = T(d) * T(0.5)
    y = abs(x) > _over_sinhcosh(T) ? typemax(T) : y
    y = isnan(y) ? typemax(T) : y
    y = isnan(x) ? T(NaN)     : y
    return y
end

# overflow argument value
_over_tanh(::Type{Float64}) = 18.714973875
_over_tanh(::Type{Float32}) = 8.664339742f0

function tanh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    e = ddrec(d)
    d = dddiv(ddsub(d, e), ddadd(d, e))
    y = T(d)
    y = abs(x) > _over_tanh(T) ? T(1) : y
    y = isnan(y) ? T(1) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end

function asinh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = logk2(ddadd(ddsqrt(ddadd2(ddsqu(y),  T(1))), y))
    y = T(d)
    y = isinf(x) || isnan(y) ? typemax(T) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end

function acosh{T<:FloatTypes}(x::T)
    d = logk2(ddadd2(ddsqrt(ddadd2(ddsqu(x), -T(1))), x))
    y = T(d)
    y = isinf(x) || isnan(y) ? typemax(T) : y
    y = x == T(1) ? T(0) : y
    y = x < T(1)  ? T(NaN)  : y
    y = isnan(x)  ? T(NaN)  : y
    return y
end

function atanh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = logk2(dddiv(ddadd2(T(1), y), ddadd2(T(1), -y)))
    y = y > T(1) ? T(NaN) : (y == T(1) ? typemax(T) : T(d) * T(0.5))
    y = isinf(x) || isnan(y) ? T(NaN) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end
