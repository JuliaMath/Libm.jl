function xsinh(x::Float64)
    y = abs(x)
    d = expk2(Double(y))
    d = ddsub(d, ddrec(d))
    y = (d.hi + d.lo) * 0.5
    y = abs(x) > 710 ? Inf : y
    y = isnan(y) ? Inf : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xcosh(x::Float64)
    y = abs(x)
    d = expk2(Double(y))
    d = ddadd(d, ddrec(d))
    y = (d.hi + d.lo) * 0.5
    y = abs(x) > 710 ? Inf : y
    y = isnan(y) ? Inf : y
    y = isnan(x) ? NaN : y
    return y
end

function xtanh(x::Float64)
    y = abs(x)
    d = expk2(Double(y))
    e = ddrec(d)
    d = dddiv(ddsub(d, e), ddadd(d, e))
    y = d.hi + d.lo
    y = abs(x) > 18.714973875 ? 1.0 : y
    y = isnan(y) ? 1.0 : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xasinh(x::Float64)
    y = abs(x)
    d = logk2(ddadd(ddsqrt(ddadd2(ddsqu(y),  1.0)), y))
    y = d.hi + d.lo
    y = isinf(x) || isnan(y) ? Inf : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xacosh(x::Float64)
    d = logk2(ddadd2(ddsqrt(ddadd2(ddsqu(x), -1.0)), x))
    y = d.hi + d.lo
    y = isinf(x) || isnan(y) ? Inf : y
    y = x == 1.0 ? 0.0 : y
    y = x < 1.0  ? NaN : y
    y = isnan(x) ? NaN : y
    return y
end

function xatanh(x::Float64)
    y = abs(x)
    d = logk2(dddiv(ddadd2(1.0, y), ddadd2(1.0, -y)))
    y = y > 1.0 ? NaN : (y == 1.0 ? Inf : (d.hi + d.lo) * 0.5)
    y = isinf(x) || isnan(y) ? NaN : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end
