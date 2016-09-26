## TODO: Clean up methods / constans / add macros

let
global xsin
global xcos

const c9d = -7.972559550090378688919520e-18
const c8d =  2.810099727108632000912510e-15
const c7d = -7.647122191181588332884840e-13
const c6d =  1.605904306056645016290540e-10
const c5d = -2.505210837635020458107550e-08
const c4d =  2.755731922391987476304160e-06
const c3d = -0.0001984126984126961628068090
const c2d =  0.0083333333333333297482381500
const c1d = -0.1666666666666666574148080000

# c5f is 0f0 to handle Inf32 case, Float64 doesn't need this since it comes
# out properly (add another neg constant and remove this zero constant)
const c5f =  0f0
const c4f =  2.608315980978659354150300f-06
const c3f = -0.00019810690719168633222580f0
const c2f =  0.00833307858556509017944336f0
const c1f = -0.16666659712791442871093800f0

# Argument is first reduced to the domain 0 < s < π/4

# We return the correct sign using `q & 1 != 0` i.e. q is odd (this works for
# positive and negative q) and if this condition is true we flip the sign since
# we are now in the negative branch of sin(x). Recall that q is just the integer
# part of d/π and thus we can determine the correct sign using this information.

global @inline _sincos(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d
global @inline _sincos(x::Float32) = @horner x c1f c2f c3f c4f c5f

function xsin{T<:FloatTypes}(x::T)
    d = abs(x)
    q = xrint(d*T(M1PI))
    d = muladd(q, -PI4A(T)*4, d)
    d = muladd(q, -PI4B(T)*4, d)
    d = muladd(q, -PI4C(T)*4, d)
    d = muladd(q, -PI4D(T)*4, d)
    s = d*d
    q & 1 != 0 && (d = -d)
    u =_sincos(s)
    u = muladd(s, u*d, d)
    return flipsign(u,x)
end

function xcos{T<:FloatTypes}(d::T)
    q = muladd(2, xrint(d*T(M1PI)-T(0.5)), 1)
    d = muladd(q, -PI4A(T)*2, d)
    d = muladd(q, -PI4B(T)*2, d)
    d = muladd(q, -PI4C(T)*2, d)
    d = muladd(q, -PI4D(T)*2, d)
    s = d*d
    q & 2 == 0 && (d = -d)
    u =_sincos(s)
    return muladd(s, u*d, d)
end
end

let
global xsin_u1
global xcos_u1

const c7d =  2.72052416138529567917983e-15
const c6d = -7.64292594113954471900203e-13
const c5d =  1.60589370117277896211623e-10
const c4d = -2.50521068148431233593680e-08
const c3d =  2.75573192104428224777379e-06
const c2d = -0.000198412698412046454654947
const c1d =  0.008333333333333180562019220

const c3f =  2.608315980978659354150300f-06
const c2f = -0.00019810690719168633222580f0
const c1f =  0.00833307858556509017944336f0

global @inline _sincos_u1(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d
global @inline _sincos_u1(x::Float32) = @horner x c1f c2f c3f

global _sincos_u1_c0(::Type{Float64}) = -0.166666666666666657414808
global _sincos_u1_c0(::Type{Float32}) = -0.166666597127914428710938f0

function xsin_u1{T<:FloatTypes}(x::T)
    d = abs(x)
    q = xrint(d*T(M1PI))
    s = ddadd2(d, q * (-PI4A(T)*4))
    s = ddadd2(s, q * (-PI4B(T)*4))
    s = ddadd2(s, q * (-PI4C(T)*4))
    s = ddadd2(s, q * (-PI4D(T)*4))
    t = s
    s = ddsqu(s)
    u =_sincos_u1(s.hi)
    v = ddadd(T(1), ddmul(ddadd(_sincos_u1_c0(T), u*s.hi), s))
    v = ddmul(t, v)
    u = v.hi + v.lo
    q & 1 != 0 && (u = -u)
    return flipsign(u,x)
end

function xcos_u1{T<:FloatTypes}(d::T)
    d = abs(d)
    q = muladd(2, xrint(d*T(M1PI) - T(0.5)), 1)
    s = ddadd2(d, q * (-PI4A(T)*2))
    s = ddadd2(s, q * (-PI4B(T)*2))
    s = ddadd2(s, q * (-PI4C(T)*2))
    s = ddadd2(s, q * (-PI4D(T)*2))
    t = s
    s = ddsqu(s)
    u =_sincos_u1(s.hi)
    v = ddadd(T(1), ddmul(ddadd(_sincos_u1_c0(T), u*s.hi), s))
    v = ddmul(t, v)
    u = v.hi + v.lo
    q & 2 == 0 && (u = -u)
    return u
end
end

let
global xsincos
global xsincos_u1

const a6d =  1.58938307283228937328511e-10
const a5d = -2.50506943502539773349318e-08
const a4d =  2.75573131776846360512547e-06
const a3d = -0.000198412698278911770864914
const a2d =  0.008333333333319184596174600
const a1d = -0.166666666666666130709393000

const a3f = -0.000195169282960705459117889f0
const a2f =  0.008332157507538795471191410f0
const a1f = -0.166666537523269653320312000f0

const b7d = -1.13615350239097429531523e-11
const b6d =  2.08757471207040055479366e-09
const b5d = -2.75573144028847567498567e-07
const b4d =  2.48015872890001867311915e-05
const b3d = -0.001388888888887140192823290
const b2d =  0.041666666666666551959206200
const b1d = -0.500000000000000000000000000

const b5f = -2.7181184236724220681935500f-07
const b4f =  2.4799044695100747048854800f-05
const b3f = -0.001388887874782085418701170f0
const b2f =  0.041666664183139801025390600f0
const b1f = -0.500000000000000000000000000f0

global @inline _sincos_a(x::Float64) = @horner x a1d a2d a3d a4d a5d a6d
global @inline _sincos_a(x::Float32) = @horner x a1f a2f a3f
global @inline _sincos_b(x::Float64) = @horner x b1d b2d b3d b4d b5d b6d b7d
global @inline _sincos_b(x::Float32) = @horner x b1f b2f b3f b4f b5f

function xsincos{T<:FloatTypes}(x::T)
    d = abs(x)
    q = xrint(d*T(M2PI))
    s = d
    s = muladd(q, -PI4A(T)*2, s)
    s = muladd(q, -PI4B(T)*2, s)
    s = muladd(q, -PI4C(T)*2, s)
    s = muladd(q, -PI4D(T)*2, s)
    t = s
    s = s*s
    u =_sincos_a(s)
    u = u * s * t
    rx = t + u
    u =_sincos_b(s)
    ry = u * s + T(1)
    q & 1 != 0     && (s = ry; ry = rx; rx = s)
    q & 2 != 0     && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d)       && (rx = ry = T(NaN))
    return Double(flipsign(rx,x), ry)
end

function xsincos_u1{T<:FloatTypes}(x::T)
    d = abs(x)
    q = xrint(d*2*T(M1PI))
    s = ddadd2(d, q * (-PI4A(T)*2))
    s = ddadd2(s, q * (-PI4B(T)*2))
    s = ddadd2(s, q * (-PI4C(T)*2))
    s = ddadd2(s, q * (-PI4D(T)*2))
    t = s
    s = ddsqu(s)
    sx = s.hi + s.lo
    u = _sincos_a(sx)
    u *= sx * t.hi
    v = ddadd(t, u)
    rx = v.hi + v.lo
    u = _sincos_b(sx)
    v = ddadd(T(1), ddmul(sx, u))
    ry = v.hi + v.lo
    q & 1 != 0     && (u = ry; ry = rx; rx = u)
    q & 2 != 0     && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d)       && (rx = ry = T(NaN))
    return Double(flipsign(rx,x),ry)
end
end

let
global xtan
global xtan_u1

const c15d =  1.01419718511083373224408e-05
const c14d = -2.59519791585924697698614e-05
const c13d =  5.23388081915899855325186e-05
const c12d = -3.05033014433946488225616e-05
const c11d =  7.14707504084242744267497e-05
const c10d =  8.09674518280159187045078e-05
const c9d  =  0.000244884931879331847054404
const c8d  =  0.000588505168743587154904506
const c7d  =  0.001456127889228124279788480
const c6d  =  0.003592087438369066191429240
const c5d  =  0.008863239443624016181133560
const c4d  =  0.021869488285384638959207800
const c3d  =  0.053968253978129841763600200
const c2d  =  0.133333333333125941821962000
const c1d  =  0.333333333333334980164153000

const c7f =  0.00446636462584137916564941f0
const c6f = -8.3920182078145444393158f-05
const c5f =  0.0109639242291450500488281f0
const c4f =  0.0212360303848981857299805f0
const c3f =  0.0540687143802642822265625f0
const c2f =  0.133325666189193725585938f0
const c1f =  0.33333361148834228515625f0

global @inline _tan(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d
global @inline _tan(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f

function xtan{T<:FloatTypes}(d::T)
    q = xrint(d*T(M2PI))
    x = muladd(q, -PI4A(T)*2, d)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    q & 1 != 0 && (x = -x)
    s = x*x
    u =_tan(s)
    u = muladd(s, u * x, x)
    q & 1 != 0 && (u = 1/u)
    isinf(d)   && (u = T(NaN))
    return u
end

global @inline _tan_u1(x::Double{Float64}) = ddadd(c1d, x.hi*(@horner x.hi c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d))
global @inline _tan_u1(x::Double{Float32}) = ddadd(c1f, ddmul(x.hi, ddadd(c2f, x.hi*(@horner x.hi c3f c4f c5f c6f c7f))))
# global @inline _tan_u1(x::Double{Float32}) = ddadd(c1f, ddmul(x, @horner x.hi c2f c3f c4f c5f c6f c7f))

function xtan_u1{T<:FloatTypes}(d::T)
    q = xrint(d*T(M2PI))
    x = ddadd2(d, q * -PI4A(T)*2)
    x = ddadd2(x, q * -PI4B(T)*2)
    x = ddadd2(x, q * -PI4C(T)*2)
    x = ddadd2(x, q * -PI4D(T)*2)
    q & 1 != 0 && (x = -x)
    s = ddsqu(x)
    u =_tan_u1(s)
    u = ddmul(x, ddadd(T(1), ddmul(u, s)))
    q & 1 != 0 && (u = ddrec(u))
    return u.hi + u.lo
end
end

let
global xatan
const c19d = -1.88796008463073496563746e-05
const c18d =  0.000209850076645816976906797
const c17d = -0.001106118314866724825634710
const c16d =  0.003700267441887131192324030
const c15d = -0.008898961958876554917408090
const c14d =  0.016599329773529201970117000
const c13d = -0.025451762493231264161686100
const c12d =  0.033785258000135306999389700
const c11d = -0.040762919127683650000193400
const c10d =  0.046666715007784062563267500
const c9d  = -0.052367485230348245761611300
const c8d  =  0.058766639292667358085431300
const c7d  = -0.066657357936108052598456200
const c6d  =  0.076921953831176961835502900
const c5d  = -0.090908995008245008229153000
const c4d  =  0.111111105648261418443745000
const c3d  = -0.142857142667713293837650000
const c2d  =  0.199999999996591265594148000
const c1d  = -0.333333333333311110369124000

const c8f =  0.00282363896258175373077393f0
const c7f = -0.01595690287649631500244140f0
const c6f =  0.04250498861074447631835940f0
const c5f = -0.07489009201526641845703120f0
const c4f =  0.10634793341159820556640600f0
const c3f = -0.14202736318111419677734400f0
const c2f =  0.19992695748805999755859400f0
const c1f = -0.33333101868629455566406200f0

global @inline _atan(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d
global @inline _atan(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f

function xatan{T<:FloatTypes}(s::T)
    q = 0
    if s < 0
        s = -s
        q = 2
    end
    if s > 1
        s = 1.0/s
        q |= 1
    end
    t = s*s
    u =_atan(t)
    t = s + s*t*u
    q & 1 != 0 && (t = T(MPI2) - t)
    q & 2 != 0 && (t = -t)
    return t
end
end

function xatan_u1{T<:FloatTypes}(d::T)
    d2 = atan2k_u1(Double(abs(d)), Double(T(1)))
    r = d2.hi + d2.lo
    isinf(d) && (r = T(MPI2))
    return flipsign(r, d)
end

function xatan2{T<:FloatTypes}(y::T, x::T)
    r = atan2k(abs(y), x)
    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI2) : T(0))
    end
    if isinf(y)
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI4) : T(0))
    end
    if y == 0
        r = _sign(x) == -1 ? T(MPI) : T(0)
    end
    return isnan(x) || isnan(y) ? T(NaN) : flipsign(r, y)
end

function xatan2_u1{T<:FloatTypes}(y::T, x::T)
    d = atan2k_u1(Double(abs(y)), Double(x))
    r = d.hi + d.lo
    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI2) : T(0))
    end
    if isinf(y)
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI4) : T(0))
    end
    if y == 0
        r = _sign(x) == -1 ? T(MPI) : T(0)
    end
    return isnan(x) || isnan(y) ? T(NaN) : flipsign(r, y)
end

xasin{T<:FloatTypes}(x::T) = flipsign(atan2k(abs(x), _sqrt((1+x)*(1-x))), x)

function xasin_u1{T<:FloatTypes}(d::T)
    d2 = atan2k_u1(Double(abs(d)), ddsqrt(ddmul(ddadd(T(1), d), ddadd(T(1),-d))))
    r = d2.hi + d2.lo
    abs(d) == 1 && (r = T(MPI2))
    return flipsign(r, d)
end

xacos{T<:FloatTypes}(x::T) = flipsign(atan2k(_sqrt((1+x)*(1-x)), abs(x)), x) + (x < 0 ? T(MPI) : T(0))

function xacos_u1{T<:FloatTypes}(d::T)
    d2 = atan2k_u1(ddsqrt(ddmul(ddadd(T(1), d), ddadd(T(1),-d))), Double(abs(d)))
    d2 = ddscale(d2, _sign(d))
    abs(d) == 1 && (d2 = Double(T(0)))
    d < 0       && (d2 = ddadd(MDPI(T), d2))
    return d2.hi + d2.lo
end
