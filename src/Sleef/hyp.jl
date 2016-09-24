# overflow argument value
_max_val_sinhcosh(::Type{Float64}) = 710.0
_max_val_sinhcosh(::Type{Float32}) = 89f0

function xsinh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    d = ddsub(d, ddrec(d))
    y = (d.hi + d.lo) * T(0.5)
    y = abs(x) > _max_val_sinhcosh(T) ? T(Inf) : y
    y = isnan(y) ? T(Inf) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end

function xcosh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    d = ddadd(d, ddrec(d))
    y = (d.hi + d.lo) * T(0.5)
    y = abs(x) > _max_val_sinhcosh(T) ? T(Inf) : y
    y = isnan(y) ? T(Inf) : y
    y = isnan(x) ? T(NaN) : y
    return y
end

# overflow argument value
_max_val_tanh(::Type{Float64}) = 18.714973875
_max_val_tanh(::Type{Float32}) = 8.664339742f0

function xtanh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = expk2(Double(y))
    e = ddrec(d)
    d = dddiv(ddsub(d, e), ddadd(d, e))
    y = d.hi + d.lo
    y = abs(x) > _max_val_tanh(T) ? T(1.0) : y
    y = isnan(y) ? T(1.0) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end

function xasinh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = logk2(ddadd(ddsqrt(ddadd2(ddsqu(y),  T(1.0))), y))
    y = d.hi + d.lo
    y = isinf(x) || isnan(y) ? T(Inf) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end

function xacosh{T<:FloatTypes}(x::T)
    d = logk2(ddadd2(ddsqrt(ddadd2(ddsqu(x), -T(1.0))), x))
    y = d.hi + d.lo
    y = isinf(x) || isnan(y) ? T(Inf) : y
    y = x == T(1.0) ? T(0.0) : y
    y = x < T(1.0)  ? T(NaN) : y
    y = isnan(x)    ? T(NaN) : y
    return y
end

function xatanh{T<:FloatTypes}(x::T)
    y = abs(x)
    d = logk2(dddiv(ddadd2(T(1.0), y), ddadd2(T(1.0), -y)))
    y = y > T(1.0) ? T(NaN) : (y == T(1.0) ? T(Inf) : (d.hi + d.lo) * T(0.5))
    y = isinf(x) || isnan(y) ? T(NaN) : y
    y = flipsign(y, x)
    y = isnan(x) ? T(NaN) : y
    return y
end
