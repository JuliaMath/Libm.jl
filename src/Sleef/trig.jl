let
global xsin
global xcos

const c9d = -7.97255955009037868891952e-18
const c8d =  2.81009972710863200091251e-15
const c7d = -7.64712219118158833288484e-13
const c6d =  1.60590430605664501629054e-10
const c5d = -2.50521083763502045810755e-08
const c4d =  2.75573192239198747630416e-06
const c3d = -0.000198412698412696162806809
const c2d =  0.00833333333333332974823815
const c1d = -0.166666666666666657414808

const c4f =  2.6083159809786593541503f-06
const c3f = -0.0001981069071916863322258f0
const c2f =  0.00833307858556509017944336f0
const c1f = -0.166666597127914428710938f0

# Argument is first reduced to the domain 0 < s < π/4

"""
We return the correct sign using `q & 1 != 0` i.e. q is odd (this works for
positive and negative q) and if this condition is true we flip the sign since
we are now in the negative branch of sin(x). Recall that q is just the integer
part of d/π and thus we can determine the correct sign using this information.
"""
global @inline _sincos(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d
global @inline _sincos(x::Float32) = @horner x c1f c2f c3f c4f

function xsin{T<:FloatTypes}(x::T)
    d = abs(x)
    q = xrint(d*M1PI)
    d = muladd(q, -PI4A(T)*4, d)
    d = muladd(q, -PI4B(T)*4, d)
    d = muladd(q, -PI4C(T)*4, d)
    d = muladd(q, -PI4D(T)*4, d)
    s = d*d
    (q & 1) != 0 && (d = -d)
    u =_sincos(s)
    u = muladd(s, u*d, d)
    return flipsign(u,x)
end

function xcos{T<:FloatTypes}(d::T)
    q = muladd(2, xrint(d*M1PI - 0.5), 1)
    d = muladd(q, -PI4A(T)*2, d)
    d = muladd(q, -PI4B(T)*2, d)
    d = muladd(q, -PI4C(T)*2, d)
    d = muladd(q, -PI4D(T)*2, d)
    s = d*d
    (q & 2) == 0 && (d = -d)
    u =_sincos(s)
    return muladd(s, u*d, d)
end
end

let
global xsin_u1
global xcos_u1

const c7 =  2.72052416138529567917983e-15
const c6 = -7.6429259411395447190023e-13
const c5 =  1.60589370117277896211623e-10
const c4 = -2.5052106814843123359368e-08
const c3 =  2.75573192104428224777379e-06
const c2 = -0.000198412698412046454654947
const c1 =  0.00833333333333318056201922

function xsin_u1(x::Float64)
    d = abs(x)
    q = xrint(d*M1PI)
    s = ddadd2(d, q * (-PI4Ad*4))
    s = ddadd2(s, q * (-PI4Bd*4))
    s = ddadd2(s, q * (-PI4Cd*4))
    s = ddadd2(s, q * (-PI4Dd*4))
    t = s
    s = ddsqu(s)
    u = @horner s.hi c1 c2 c3 c4 c5 c6 c7
    v = ddadd(1.0, ddmul(ddadd(-0.166666666666666657414808, u*s.hi), s))
    v = ddmul(t, v)
    u = v.hi + v.lo
    (q & 1) != 0 && (u = -u)
    return flipsign(u,x)
end

function xcos_u1(d::Float64)
    d = abs(d)
    q = muladd(2, xrint(d * M1PI - 0.5), 1)
    s = ddadd2(d, q * (-PI4Ad*2))
    s = ddadd2(s, q * (-PI4Bd*2))
    s = ddadd2(s, q * (-PI4Cd*2))
    s = ddadd2(s, q * (-PI4Dd*2))
    t = s
    s = ddsqu(s)
    u = @horner s.hi c1 c2 c3 c4 c5 c6 c7
    v = ddadd(1.0, ddmul(ddadd(-0.166666666666666657414808, u*s.hi), s))
    v = ddmul(t, v)
    u = v.hi + v.lo
    (q & 2) == 0 && (u = -u)
    return u
end
end

let
global xsincos
global xsincos_u1

const a6 =  1.58938307283228937328511e-10
const a5 = -2.50506943502539773349318e-08
const a4 =  2.75573131776846360512547e-06
const a3 = -0.000198412698278911770864914
const a2 =  0.0083333333333191845961746
const a1 = -0.166666666666666130709393

const b7 = -1.13615350239097429531523e-11
const b6 =  2.08757471207040055479366e-09
const b5 = -2.75573144028847567498567e-07
const b4 =  2.48015872890001867311915e-05
const b3 = -0.00138888888888714019282329
const b2 =  0.0416666666666665519592062
const b1 = -0.5

function xsincos(x::Float64)
    d = abs(x)
    q = xrint(d*(2*M1PI))
    s = d
    s = muladd(-q, PI4Ad*2, s)
    s = muladd(-q, PI4Bd*2, s)
    s = muladd(-q, PI4Cd*2, s)
    s = muladd(-q, PI4Dd*2, s)
    t = s
    s = s*s
    u = @horner s a1 a2 a3 a4 a5 a6
    u = u * s * t
    rx = t + u
    u = @horner s b1 b2 b3 b4 b5 b6 b7
    ry = u * s + 1
    (q & 1) != 0     && (s = ry; ry = rx; rx = s)
    (q & 2) != 0     && (rx = -rx)
    ((q+1) & 2) != 0 && (ry = -ry)
    isinf(d) && (rx = ry = NaN)
    return Double(flipsign(rx,x),ry)
end

function xsincos_u1(x::Float64)
    d = abs(x)
    q = xrint(d*(2*M1PI))
    s = ddadd2(d, q * (-PI4Ad*2))
    s = ddadd2(s, q * (-PI4Bd*2))
    s = ddadd2(s, q * (-PI4Cd*2))
    s = ddadd2(s, q * (-PI4Dd*2))
    t = s
    s = ddsqu(s)
    sx = s.hi + s.lo
    u = @horner sx a1 a2 a3 a4 a5 a6
    u *= sx * t.hi
    v = ddadd(t, u)
    rx = v.hi + v.lo
    u = @horner sx b1 b2 b3 b4 b5 b6 b7
    v = ddadd(1.0, ddmul(sx, u))
    ry = v.hi + v.lo
    (q & 1) != 0     && (u = ry; ry = rx; rx = u)
    (q & 2) != 0     && (rx = -rx)
    ((q+1) & 2) != 0 && (ry = -ry)
    isinf(d) && (rx = ry = NaN)
    return Double(flipsign(rx,x),ry)
end
end

let
global xtan
global xtan_u1

const c15 =  1.01419718511083373224408e-05
const c14 = -2.59519791585924697698614e-05
const c13 =  5.23388081915899855325186e-05
const c12 = -3.05033014433946488225616e-05
const c11 =  7.14707504084242744267497e-05
const c10 =  8.09674518280159187045078e-05
const c9  =  0.000244884931879331847054404
const c8  =  0.000588505168743587154904506
const c7  =  0.00145612788922812427978848
const c6  =  0.00359208743836906619142924
const c5  =  0.00886323944362401618113356
const c4  =  0.0218694882853846389592078
const c3  =  0.0539682539781298417636002
const c2  =  0.133333333333125941821962
const c1  =  0.333333333333334980164153

function xtan(d::Float64)
    q = xrint(d * (2*M1PI))
    x = muladd(q, -PI4Ad*2, d)
    x = muladd(q, -PI4Bd*2, x)
    x = muladd(q, -PI4Cd*2, x)
    x = muladd(q, -PI4Dd*2, x)
    s = x*x
    (q & 1) != 0 && (x = -x)
    u = @horner s c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15
    u = muladd(s, u * x, x)
    (q & 1) != 0 && (u = 1.0/u)
    isinf(d) && (u = NaN)
    return u
end

function xtan_u1(d::Float64)
    q = xrint(d*M2PI)
    s = ddadd2(d, q*(-PI4Ad*2))
    s = ddadd2(s, q*(-PI4Bd*2))
    s = ddadd2(s, q*(-PI4Cd*2))
    s = ddadd2(s, q*(-PI4Dd*2))
    (q & 1) != 0 && (s = ddneg(s))
    t = s
    s = ddsqu(s)
    u = @horner s.hi c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15
    x = ddadd(1.0, ddmul(ddadd(0.333333333333334980164153, u * s.hi), s))
    x = ddmul(t, x)
    (q & 1) != 0 && (x = ddrec(x))
    u = x.hi + x.lo
    return u
end
end

function xatan2(y::Float64, x::Float64)
    r = atan2k(abs(y), x)
    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = MPI/2 - (isinf(x) ? (_sign(x) * (MPI/2)) : 0.0)
    end
    if isinf(y)
        r = MPI/2 - (isinf(x) ? (_sign(x) * (MPI/4)) : 0.0)
    end
    if y == 0
        r = (_sign(x) == -1 ? MPI/1 : 0.0)
    end
    return isnan(x) || isnan(y) ? NaN : flipsign(r, y)
end

xasin(d::Float64) = flipsign(atan2k(abs(d), _sqrt((1+d)*(1-d))), d)

xacos(d::Float64) = flipsign(atan2k(_sqrt((1+d)*(1-d)), abs(d)), d) + (d < 0 ? MPI/1 : 0.0)

let # check constants and figure out why @horner fails here
global xatan
const c19 = -1.88796008463073496563746e-05
const c18 =  0.000209850076645816976906797
const c17 = -0.00110611831486672482563471
const c16 =  0.00370026744188713119232403
const c15 = -0.00889896195887655491740809
const c14 =  0.016599329773529201970117
const c13 = -0.0254517624932312641616861
const c12 =  0.0337852580001353069993897
const c11 = -0.0407629191276836500001934
const c10 =  0.0466667150077840625632675
const c9  = -0.0523674852303482457616113
const c8  =  0.0587666392926673580854313
const c7  = -0.0666573579361080525984562
const c6  =  0.0769219538311769618355029
const c5  = -0.090908995008245008229153
const c4  =  0.111111105648261418443745
const c3  = -0.14285714266771329383765
const c2  =  0.199999999996591265594148
const c1  = -0.333333333333311110369124

function xatan(s::Float64)
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
    u = c19
    u = u * t + c18
    u = u * t + c17
    u = u * t + c16
    u = u * t + c15
    u = u * t + c14
    u = u * t + c13
    u = u * t + c12
    u = u * t + c11
    u = u * t + c10
    u = u * t + c9
    u = u * t + c8
    u = u * t + c7
    u = u * t + c6
    u = u * t + c5
    u = u * t + c4
    u = u * t + c3
    u = u * t + c2
    u = u * t + c1
    t = s + s*(t*u)
    (q & 1) != 0 && (t = 1.570796326794896557998982 - t)
    (q & 2) != 0 && (t = -t)
    return t
end
end

function xatan2_u1(y::Float64, x::Float64)
    d = atan2k_u1(Double(abs(y)), Double(x))
    r = d.hi + d.lo

    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = MPI/2 - (isinf(x) ? (_sign(x) * (MPI/2)) : 0.0)
    end
    if isinf(y)
        r = MPI/2 - (isinf(x) ? (_sign(x) * (MPI/4)) : 0.0)
    end
    if y == 0
        r = _sign(x) == -1 ? MPI/1 : 0.0
    end
    return isnan(x) || isnan(y) ? NaN : flipsign(r, y)
end

function xasin_u1(d::Float64)
    d2 = atan2k_u1(Double(abs(d)), ddsqrt(ddmul(ddadd(1.0, d), ddadd(1.0,-d))))
    r = d2.hi + d2.lo
    abs(d) == 1 && (r = 1.570796326794896557998982)
    return flipsign(r, d)
end

function xacos_u1(d::Float64)
    d2 = atan2k_u1(ddsqrt(ddmul(ddadd(1.0, d), ddadd(1.0,-d))), Double(abs(d)))
    d2 = ddscale(d2, _sign(d))
    abs(d) == 1 && (d2 = Double(0.0))
    d < 0 && (d2 = ddadd(Double(3.141592653589793116, 1.2246467991473532072e-16), d2))
    return d2.hi + d2.lo
end

function xatan_u1(d::Float64)
    d2 = atan2k_u1(Double(abs(d)), Double(1.0))
    r = d2.hi + d2.lo
    isinf(d) && (r = 1.570796326794896557998982)
    return flipsign(r, d)
end
