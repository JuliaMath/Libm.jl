
function xpow(x::Float64, y::Float64)
    yint = unsafe_trunc(Int32,y)
    yisint = yint == y
    yisodd = isodd(yint) && yisint

    result = expk(ddmul_d2_d2_d(logk(abs(x)), y))

    result = isnan(result) ? Inf : result
    result *=  (x >= 0 ? 1.0 : (!yisint ? NaN : (yisodd ? -1 : 1)));

    efx = flipsign(abs(x) - 1, y)
    if isinf(y)
        result = efx < 0 ? 0.0 : (efx == 0 ? 1.0 : Inf)
    end
    if isinf(x) || x == 0
        result = (yisodd ? sign(x) : 1.0) * ((x == 0 ? -y : y) < 0 ? 0.0 : Inf)
    end
    (isnan(x) || isnan(y)) && (result = NaN)
    (y == 0 || x == 1) && (result = 1.0)
    return result
end


function xcbrt(d::Float64) # max error 2 ulps
    q = 1.0

    e = ilogbp1(d)
    d = ldexpk(d, -e)
    r = (e + 6144) % 3
    q = (r == 1) ? 1.2599210498948731647672106 : q
    q = (r == 2) ? 1.5874010519681994747517056 : q
    q = ldexpk(q, (e + 6144)÷3 - 2048)

    q = flipsign(q, d)
    d = abs(d)

    x = -0.640245898480692909870982
    x = x * d + 2.96155103020039511818595
    x = x * d + -5.73353060922947843636166
    x = x * d + 6.03990368989458747961407
    x = x * d + -3.85841935510444988821632
    x = x * d + 2.2307275302496609725722

    y = x*x
    y = y*y
    x -= (d*y - x)*(1.0/3.0)
    y = d*x*x
    y = (y - (2.0/3.0)*y*(y*x - 1))*q

    return y
end

function xcbrt_u1(d::Float64)
    q2 = Double(1.0, 0.0)
    
    e = ilogbp1(d)
    d = ldexpk(d, -e)
    r = (e + 6144) % 3
    q2 = (r == 1) ? Double(1.2599210498948731907, -2.5899333753005069177e-17) : q2
    q2 = (r == 2) ? Double(1.5874010519681995834, -1.0869008194197822986e-16) : q2

    q3 = Double(flipsign(q2.x, d), flipsign(q2.y, d))
    d = abs(d)
    
    x = -0.640245898480692909870982
    x = x * d + 2.96155103020039511818595
    x = x * d + -5.73353060922947843636166
    x = x * d + 6.03990368989458747961407
    x = x * d + -3.85841935510444988821632
    x = x * d + 2.2307275302496609725722

    y = x*x
    y = y*y
    x -= (d*y - x)*(1.0/3.0)

    z = x

    u = ddmul_d2_d_d(x, x)
    u = ddmul_d2_d2_d2(u, u)
    u = ddmul_d2_d2_d(u, d)
    u = ddadd2_d2_d2_d(u, -x)
    y = u.x + u.y

    y = -2.0/3.0*y*z
    v = ddadd2_d2_d2_d(ddmul_d2_d_d(z, z), y)
    v = ddmul_d2_d2_d(v, d)
    v = ddmul_d2_d2_d2(v, q3)
    z = ldexp(v.x + v.y, (e + 6144)÷3 - 2048)

    isinf(d) && (z = flipsign(Inf, q3.x))
    d == 0 && (z = flipsign(0.0, q3.x))

    return z
end