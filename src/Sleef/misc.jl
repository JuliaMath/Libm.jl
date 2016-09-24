
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
const c6d = -0.640245898480692909870982
const c5d = 2.96155103020039511818595
const c4d = -5.73353060922947843636166
const c3d = 6.03990368989458747961407
const c2d = -3.85841935510444988821632
const c1d = 2.2307275302496609725722

const c6f = -0.601564466953277587890625f0
const c5f =  2.8208892345428466796875f0
const c4f = -5.532182216644287109375f0
const c3f =  5.898262500762939453125f0
const c2f = -3.8095417022705078125f0
const c1f =  2.2241256237030029296875f0

global @inline _xcbrt(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d
global @inline _xcbrt(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f

global _xcbrt_u1_a0(::Type{Float64}) = Double(1.2599210498948731907, -2.5899333753005069177e-17)
global _xcbrt_u1_a1(::Type{Float64}) = Double(1.5874010519681995834, -1.0869008194197822986e-16)
global _xcbrt_u1_a0(::Type{Float32}) = Double(1.2599210739135742188f0, -2.4018701694217270415f-08)
global _xcbrt_u1_a1(::Type{Float32}) = Double(1.5874010324478149414f0,  1.9520385308169352356f-08)


function xcbrt{T<:FloatTypes}(d::T) # max error 2 ulps
    q  = T(1)
    e  = ilogbp1(d)
    d  = ldexpk(d, -e)
    r  = (e + 6144) % 3
    q  = (r == 1) ? T(1.2599210498948731647672106) : q
    q  = (r == 2) ? T(1.5874010519681994747517056) : q
    q  = ldexpk(q, (e + 6144)÷3 - 2048)
    q  = flipsign(q, d)
    d  = abs(d)
    x  =_xcbrt(d)
    y  = x*x
    y  = y*y
    x -= (d*y - x)*(T(1)/T(3))
    y  = d*x*x
    y  = (y - (T(2)/T(3))*y*(y*x - 1))*q
    return y
end

function xcbrt_u1{T<:FloatTypes}(d::T)
    q2 = Double(T(1))  
    e  = ilogbp1(d)
    d  = ldexpk(d, -e)
    r  = (e + 6144) % 3
    q2 = (r == 1) ? _xcbrt_u1_a0(T) : q2
    q2 = (r == 2) ? _xcbrt_u1_a1(T) : q2
    q3 = Double(flipsign(q2.hi, d), flipsign(q2.lo, d))
    d  = abs(d)
    x  =_xcbrt(d)
    y  = x*x
    y  = y*y
    x -= (d*y - x)*(T(1)/T(3))
    z  = x
    u  = ddmul(x, x)
    u  = ddmul(u, u)
    u  = ddmul(u, d)
    u  = ddadd2(u, -x)
    y  = u.hi + u.lo
    y  = -(T(2)/T(3))*y*z
    v  = ddadd2(ddmul(z, z), y)
    v  = ddmul(v, d)
    v  = ddmul(v, q3)
    z  = ldexp(v.hi + v.lo, (e + 6144)÷3 - 2048)
    isinf(d) && (z = flipsign(T(Inf), q3.hi))
    d == 0   && (z = flipsign(T(0), q3.hi))
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