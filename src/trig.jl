

let
global sin_fast
global cos_fast

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

global @inline _sincos_fast(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d
global @inline _sincos_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f

function sin_fast{T<:FloatTypes}(x::T)
    d = abs(x)
    q = rint(d*T(M1PI))
    d = muladd(q, -PI4A(T)*4, d)
    d = muladd(q, -PI4B(T)*4, d)
    d = muladd(q, -PI4C(T)*4, d)
    d = muladd(q, -PI4D(T)*4, d)
    s = d*d
    q & 1 != 0 && (d = -d)
    u =_sincos_fast(s)
    u = muladd(s, u*d, d)
    return flipsign(u,x)
end


function cos_fast{T<:FloatTypes}(x::T)
    q = muladd(2, rint(x*T(M1PI)-T(0.5)), 1)
    x = muladd(q, -PI4A(T)*2, x)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    s = x*x
    q & 2 == 0 && (x = -x)
    u =_sincos_fast(s)
    return muladd(s, u*x, x)
end
end


let
global sin
global cos

const c8d =  2.72052416138529567917983e-15
const c7d = -7.64292594113954471900203e-13
const c6d =  1.60589370117277896211623e-10
const c5d = -2.50521068148431233593680e-08
const c4d =  2.75573192104428224777379e-06
const c3d = -0.000198412698412046454654947
const c2d =  0.008333333333333180562019220
const c1d = -0.166666666666666657414808

const c4f =  2.608315980978659354150300f-06
const c3f = -0.00019810690719168633222580f0
const c2f =  0.00833307858556509017944336f0
const c1f = -0.166666597127914428710938f0

global @inline _sincos(x::Double{Float64}) = dadd(c1d, x.hi*(@horner x.hi c2d c3d c4d c5d c6d c7d c8d))
global @inline _sincos(x::Double{Float32}) = dadd(c1f, x.hi*(@horner x.hi c2f c3f c4f))

function sin{T<:FloatTypes}(x::T)
    d = abs(x)
    q = rint(d*T(M1PI))
    s = dadd2(d, q * -PI4A(T)*4)
    s = dadd2(s, q * -PI4B(T)*4)
    s = dadd2(s, q * -PI4C(T)*4)
    s = dadd2(s, q * -PI4D(T)*4)
    t = s
    s = ddsqu(s)
    w =_sincos(s)
    v = dmul(t, dadd(T(1), dmul(w, s)))
    u = T(v)
    q & 1 != 0 && (u = -u)
    return flipsign(u,x)
end


function cos{T<:FloatTypes}(x::T)
    x = abs(x)
    q = muladd(2, rint(x*T(M1PI) - T(0.5)), 1)
    s = dadd2(x, q * -PI4A(T)*2)
    s = dadd2(s, q * -PI4B(T)*2)
    s = dadd2(s, q * -PI4C(T)*2)
    s = dadd2(s, q * -PI4D(T)*2)
    t = s
    s = ddsqu(s)
    w =_sincos(s)
    v = dmul(t, dadd(T(1), dmul(w, s)))
    u = T(v)
    q & 2 == 0 && (u = -u)
    return u
end
end


let
global sincos_fast
global sincos

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

function sincos_fast{T<:FloatTypes}(x::T)
    d  = abs(x)
    q  = rint(d*T(M2PI))
    s  = d
    s  = muladd(q, -PI4A(T)*2, s)
    s  = muladd(q, -PI4B(T)*2, s)
    s  = muladd(q, -PI4C(T)*2, s)
    s  = muladd(q, -PI4D(T)*2, s)
    t  = s
    s  = s*s
    u  =_sincos_a(s)
    u  = u * s * t
    rx = t + u
    u =_sincos_b(s)
    ry = u * s + T(1)
    q & 1 != 0     && (s = ry; ry = rx; rx = s)
    q & 2 != 0     && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d)       && (rx = ry = T(NaN))
    return Double(flipsign(rx,x), ry)
end


function sincos{T<:FloatTypes}(x::T)
    d  = abs(x)
    q  = rint(d*2*T(M1PI))
    s  = dadd2(d, q * -PI4A(T)*2)
    s  = dadd2(s, q * -PI4B(T)*2)
    s  = dadd2(s, q * -PI4C(T)*2)
    s  = dadd2(s, q * -PI4D(T)*2)
    t  = s
    s  = ddsqu(s)
    sx = T(s)
    u  =_sincos_a(sx)
    u *= sx * t.hi
    v  = dadd(t, u)
    rx = T(v)
    u  =_sincos_b(sx)
    v  = dadd(T(1), dmul(sx, u))
    ry = T(v)
    q & 1 != 0     && (u  =  ry; ry = rx; rx = u)
    q & 2 != 0     && (rx = -rx)
    (q+1) & 2 != 0 && (ry = -ry)
    isinf(d)       && (rx =  ry = T(NaN))
    return Double(flipsign(rx, x), ry)
end
end


let
global tan_fast
global tan

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

global @inline _tan_fast(x::Float64) = @horner_split x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d
global @inline _tan_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f

function tan_fast{T<:FloatTypes}(d::T)
    q = rint(d*T(M2PI))
    x = muladd(q, -PI4A(T)*2, d)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    q & 1 != 0 && (x = -x)
    s = x*x
    u =_tan_fast(s)
    u = muladd(s, u * x, x)
    q & 1 != 0 && (u = 1/u)
    isinf(d)   && (u = T(NaN))
    return u
end

global @inline _tan(x::Double{Float64}) = dadd(c1d, x.hi*(@horner_split x.hi c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d))
global @inline _tan(x::Double{Float32}) = dadd(c1f, dmul(x, @horner x.hi c2f c3f c4f c5f c6f c7f))
# global @inline _tan(x::Double{Float32}) = dadd(c1f, dmul(x.hi, dadd(c2f, x.hi*(@horner x.hi c3f c4f c5f c6f c7f))))

function tan{T<:FloatTypes}(d::T)
    q = rint(d*T(M2PI))
    x = dadd2(d, q * -PI4A(T)*2)
    x = dadd2(x, q * -PI4B(T)*2)
    x = dadd2(x, q * -PI4C(T)*2)
    x = dadd2(x, q * -PI4D(T)*2)
    q & 1 != 0 && (x = -x)
    s = ddsqu(x)
    u =_tan(s)
    u = dmul(x, dadd(T(1), dmul(u, s)))
    q & 1 != 0 && (u = ddrec(u))
    return T(u)
end
end


let
global atan_fast

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

global @inline _atan_fast(x::Float64) = @horner_split x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d
global @inline _atan_fast(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f

function atan_fast{T<:FloatTypes}(x::T)
    q = 0
    if x < 0
        x = -x
        q = 2
    end
    if x > 1
        x = 1.0/x
        q |= 1
    end
    t = x*x
    u =_atan_fast(t)
    t = x + x*t*u
    q & 1 != 0 && (t = T(MPI2) - t)
    q & 2 != 0 && (t = -t)
    return t
end
end


function atan{T<:FloatTypes}(x::T)
    x2 = atan2k(Double(abs(x)), Double(T(1)))
    r = T(x2)
    isinf(x) && (r = T(MPI2))
    return flipsign(r, x)
end


function atan2_fast{T<:FloatTypes}(y::T, x::T)
    r = atan2k_fast(abs(y), x)
    r = flipsign(r,x)
    if isinf(x) || x == 0
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI2) : T(0))
    end
    if isinf(y)
        r = T(MPI2) - (isinf(x) ? _sign(x)*T(MPI4) : T(0))
    end
    if y == 0
        r = _sign(x) == -1 ? T(MPI) : T(0)
    end
    return isnan(x) || isnan(y) ? T(NaN) : flipsign(r,y)
end


function atan2{T<:FloatTypes}(y::T, x::T)
    d = atan2k(Double(abs(y)), Double(x))
    r = T(d)
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
    return isnan(x) || isnan(y) ? T(NaN) : flipsign(r,y)
end


asin_fast{T<:FloatTypes}(x::T) = flipsign(atan2k_fast(abs(x), _sqrt((1+x)*(1-x))), x)


function asin{T<:FloatTypes}(x::T)
    x2 = atan2k(Double(abs(x)), dsqrt(dmul(dadd(T(1), x), dadd(T(1),-x))))
    u = T(x2)
    abs(x) == 1 && (u = T(MPI2))
    return flipsign(u,x)
end


acos_fast{T<:FloatTypes}(x::T) = flipsign(atan2k_fast(_sqrt((1+x)*(1-x)), abs(x)), x) + (x < 0 ? T(MPI) : T(0))


function acos{T<:FloatTypes}(x::T)
    x2 = atan2k(dsqrt(dmul(dadd(T(1), x), dsub(T(1), x))), Double(abs(x)))
    x2 = flipsign(x2,x)
    abs(x) == 1 && (x2 = Double(T(0)))
    x < 0       && (x2 = dadd(MDPI(T), x2))
    return T(x2)
end
