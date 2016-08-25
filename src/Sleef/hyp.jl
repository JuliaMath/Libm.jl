function xsinh(x::Float64)
    y = abs(x)
    d = expk2(Double2(y, 0.0))
    d = ddsub_d2_d2_d2(d, ddrec_d2_d2(d))
    y = (d.x + d.y) * 0.5

    y = abs(x) > 710 ? Inf : y
    y = isnan(y) ? Inf : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xcosh(x::Float64)
    y = abs(x)
    d = expk2(Double2(y, 0.0))
    d = ddadd_d2_d2_d2(d, ddrec_d2_d2(d))
    y = (d.x + d.y) * 0.5

    y = abs(x) > 710 ? Inf : y
    y = isnan(y) ? Inf : y
    y = isnan(x) ? NaN : y
    return y
end

function xtanh(x::Float64)
    y = abs(x)
    d = expk2(Double2(y, 0.0))
    e = ddrec_d2_d2(d)
    d = dddiv_d2_d2_d2(ddsub_d2_d2_d2(d, e), ddadd_d2_d2_d2(d, e))
    y = d.x + d.y

    y = abs(x) > 18.714973875 ? 1.0 : y
    y = isnan(y) ? 1.0 : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xasinh(x::Float64)
    y = abs(x)
    d = logk2(ddadd_d2_d2_d(ddsqrt_d2_d2(ddadd2_d2_d2_d(ddmul_d2_d_d(y, y),  1.0)), y))
    y = d.x + d.y

    y = isinf(x) || isnan(y) ? Inf : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end

function xacosh(x::Float64)
    d = logk2(ddadd2_d2_d2_d(ddsqrt_d2_d2(ddadd2_d2_d2_d(ddmul_d2_d_d(x, x), -1.0)), x))
    y = d.x + d.y

    y = isinf(x) || isnan(y) ? Inf : y
    y = x == 1.0 ? 0.0 : y
    y = x < 1.0 ? NaN : y
    y = isnan(x) ? NaN : y
    return y
end

function xatanh(x::Float64)
    y = abs(x)
    d = logk2(dddiv_d2_d2_d2(ddadd2_d2_d_d(1.0, y), ddadd2_d2_d_d(1.0, -y)))
    y = y > 1.0 ? NaN : (y == 1.0 ? Inf : (d.x + d.y) * 0.5)

    y = isinf(x) || isnan(y) ? NaN : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y
    return y
end
