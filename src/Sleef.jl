module Sleef
using Base.Math.@horner

export xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
    xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp

# higher accuracy functions
export xatan2_u1, xasin_u1, xacos_u1, xatan_u1, xsin_u1, xcos_u1, xsincos_u1, xtan_u1, xcbrt_u1, xlog_u1

# alias for supported floating point types
typealias FTypes Union{Float32,Float64}

const PI4_A = 0.78539816290140151978
const PI4_B = 4.9604678871439933374e-10
const PI4_C = 1.1258708853173288931e-18
const PI4_D = 1.7607799325916000908e-27

const M_4_PI = 1.273239544735162542821171882678754627704620361328125

const L2U = .69314718055966295651160180568695068359375
const L2L = .28235290563031577122588448175013436025525412068e-12
const R_LN2 = 1.442695040888963407359924681001892137426645954152985934135449406931

const M_1_PI = 1/pi
const M_PI = pi

include("Sleef/double2.jl")
# private math functions
include("Sleef/helper.jl")

# exported math functions
include("Sleef/exp.jl")
include("Sleef/log.jl")
include("Sleef/trig.jl")
include("Sleef/hyp.jl")

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
    q = ldexpk(q, unsafe_trunc(Int32, (e + 6144)/3 - 2048))

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
    q2 = Double2(1.0, 0.0)
    
    e = ilogbp1(d)
    d = ldexpk(d, -e)
    r = (e + 6144) % 3
    q2 = (r == 1) ? Double2(1.2599210498948731907, -2.5899333753005069177e-17) : q2
    q2 = (r == 2) ? Double2(1.5874010519681995834, -1.0869008194197822986e-16) : q2

    q3 = Double2(flipsign(q2.x, d), flipsign(q2.y, d))
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
    z = ldexp(v.x + v.y, unsafe_trunc(Int32, (e + 6144)/3 - 2048))

    isinf(d) && (z = flipsign(Inf, q3.x))
    d == 0 && (z = flipsign(0.0, q3.x))

    return z
end

end