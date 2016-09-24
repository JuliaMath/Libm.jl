
function xpow{T<:FloatTypes}(x::T, y::T)
    yint = unsafe_trunc(Int,y)
    yisint = yint == y
    yisodd = isodd(yint) && yisint

    result = expk(ddmul(logk(abs(x)), y))
    result = isnan(result) ? T(Inf) : result
    result *= (x >= 0 ? 1 : (!yisint ? T(NaN) : (yisodd ? -1 : 1)));
    efx = flipsign(abs(x) - 1, y)
    if isinf(y)
        result = efx < 0 ? T(0) : (efx == 0 ? T(1) : T(Inf))
    end
    if isinf(x) || x == 0
        result = (yisodd ? _sign(x) : T(1)) * ((x == 0 ? -y : y) < 0 ? T(0) : T(Inf))
    end
    (isnan(x) || isnan(y)) && (result = T(NaN))
    (y == 0   || x == 1)   && (result = T(1))
    return result
end

let
global xcbrt
global xcbrt_u1
const c6 = -0.640245898480692909870982
const c5 = 2.96155103020039511818595
const c4 = -5.73353060922947843636166
const c3 = 6.03990368989458747961407
const c2 = -3.85841935510444988821632
const c1 = 2.2307275302496609725722

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
    x = @horner d c1 c2 c3 c4 c5 c6
    y = x*x
    y = y*y
    x -= (d*y - x)*(1.0/3.0)
    y = d*x*x
    y = (y - (2.0/3.0)*y*(y*x - 1))*q
    return y
end

function xcbrt_u1(d::Float64)
    q2 = Double(1.0)  
    e = ilogbp1(d)
    d = ldexpk(d, -e)
    r = (e + 6144) % 3
    q2 = (r == 1) ? Double(1.2599210498948731907, -2.5899333753005069177e-17) : q2
    q2 = (r == 2) ? Double(1.5874010519681995834, -1.0869008194197822986e-16) : q2
    q3 = Double(flipsign(q2.hi, d), flipsign(q2.lo, d))
    d = abs(d)
    x = @horner d c1 c2 c3 c4 c5 c6
    y = x*x
    y = y*y
    x -= (d*y - x)*(1.0/3.0)
    z = x
    u = ddmul(x, x)
    u = ddmul(u, u)
    u = ddmul(u, d)
    u = ddadd2(u, -x)
    y = u.hi + u.lo
    y = -2.0/3.0*y*z
    v = ddadd2(ddmul(z, z), y)
    v = ddmul(v, d)
    v = ddmul(v, q3)
    z = ldexp(v.hi + v.lo, (e + 6144)÷3 - 2048)
    isinf(d) && (z = flipsign(Inf, q3.hi))
    d == 0   && (z = flipsign(0.0, q3.hi))
    return z
end
end

# custom, not in sleef, uncomment after tests added
# function xxhypot{T<:FloatTypes}(x::T, y::T)
#     x = abs(x)
#     y = abs(y)
#     if x < y
#        x, y = y, x
#     end
#     r = (x == 0) ? y : y/x
#     return x*sqrt(one(T) + r*r)
# end