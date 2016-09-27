
function xpow{T<:FloatTypes}(x::T, y::T)
    yint = unsafe_trunc(Int,y)
    yisint = yint == y
    yisodd = isodd(yint) && yisint
    z = expk(ddmul(logk(abs(x)), y))
    z = isnan(z) ? typemax(T) : z
    z *= (x >= 0 ? 1 : (!yisint ? T(NaN) : (yisodd ? -1 : 1)));
    efx = flipsign(abs(x) - 1, y)
    isinf(y) && (z = efx < 0 ? T(0) : (efx == 0 ? T(1) : typemax(T)))
    if isinf(x) || x == 0
        z = (yisodd ? _sign(x) : T(1)) * ((x == 0 ? -y : y) < 0 ? T(0) : typemax(T))
    end
    (isnan(x) || isnan(y)) && (z = T(NaN))
    (y == 0   || x == 1)   && (z = T(1))
    return z
end

let
global xcbrt_fast
global xcbrt

const c6d = -0.640245898480692909870982
const c5d = 2.9615510302003951181859500
const c4d = -5.733530609229478436361660
const c3d = 6.0399036898945874796140700
const c2d = -3.858419355104449888216320
const c1d = 2.2307275302496609725722000

const c6f = -0.601564466953277587890625f0
const c5f =  2.820889234542846679687500f0
const c4f = -5.532182216644287109375000f0
const c3f =  5.898262500762939453125000f0
const c2f = -3.809541702270507812500000f0
const c1f =  2.224125623703002929687500f0

global @inline _xcbrt(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d
global @inline _xcbrt(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f

function xcbrt_fast{T<:FloatTypes}(d::T) # max error 2 ulps
    q  = T(1)
    e  = ilogbp1(d)
    d  = ldexpk(d, -e)
    r  = (e + 6144) % 3
    q  = (r == 1) ? T(M2P13) : q
    q  = (r == 2) ? T(M2P23) : q
    q  = ldexpk(q, (e + 6144)÷3 - 2048)
    q  = flipsign(q, d)
    d  = abs(d)
    x  =_xcbrt(d)
    y  = x*x
    y  = y*y
    x -= (d*y - x)*T(1/3)
    y  = d*x*x
    y  = (y - T(2/3)*y*(y*x - 1))*q
    return y
end

function xcbrt{T<:FloatTypes}(d::T)
    q2 = Double(T(1))  
    e  = ilogbp1(d)
    d  = ldexpk(d, -e)
    r  = (e + 6144) % 3
    q2 = (r == 1) ? MD2P13(T) : q2
    q2 = (r == 2) ? MD2P23(T) : q2
    q2 = flipsign(q2, d)
    d  = abs(d)
    x  =_xcbrt(d)
    y  = x*x
    y  = y*y
    x -= (d*y - x)*T(1/3)
    z  = x
    u  = ddsqu(x)
    u  = ddsqu(u)
    u  = ddmul(u, d)
    u  = ddadd2(u, -x)
    y  = T(u)
    y  = -T(2/3)*y*z
    v  = ddadd2(ddsqu(z), y)
    v  = ddmul(v, d)
    v  = ddmul(v, q2)
    z  = ldexp(T(v), (e + 6144)÷3 - 2048)
    isinf(d) && (z = flipsign(typemax(T), q2.hi))
    d == 0   && (z = flipsign(T(0), q2.hi))
    return z
end
end

# custom, not in sleef, uncomment after tests added
# function xhypot{T<:FloatTypes}(x::T, y::T)
#     x = abs(x)
#     y = abs(y)
#     if x < y
#        x, y = y, x
#     end
#     r = (x == 0) ? y : y/x
#     return x*sqrt(T(1) + r*r)
# end