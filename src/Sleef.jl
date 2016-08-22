module Sleef
using Base.Math.@horner

export xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
    xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp

# higher accuracy functions
export xatan2_u1, xasin_u1, xacos_u1, xatan_u1, xsin_u1, xcos_u1, xsincos_u1, xtan_u1, xcbrt_u1, xlog_u1

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
include("Sleef/helper.jl")

# exported math functions

xldexp(x::Float64, q::Int32) = ldexpk(x, q)

function xilogb(d::Float64)
    e = ilogbp1(abs(d)) - 1
    e = d == 0 ? -2147483648 : e
    e = d == Inf || d == -Inf ? 2147483647 : e
    return e % Int32
end

function xatan2(y::Float64, x::Float64)
    r = atan2k(abs(y), x)

    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI /2)) : 0)
    end
    if isinf(y)
        r = M_PI/2 - (xisinf(x) ? (sign(x) * (M_PI*1/4)) : 0)
    end
    if y == 0
        r = (sign(x) == -1 ? M_PI : 0)
    end
    return isnan(x) || isnan(y) ? NaN : flipsign(r, y)
end

xasin(d::Float64) = flipsign(atan2k(abs(d), sqrt((1+d)*(1-d))), d)

xacos(d::Float64) = flipsign(atan2k(sqrt((1+d)*(1-d)), abs(d)), d) + (d < 0 ? M_PI : 0);

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

    # sleef does not use mla here
    u = -1.88796008463073496563746e-05
    u = mla(u, t, 0.000209850076645816976906797)
    u = mla(u, t, -0.00110611831486672482563471)
    u = mla(u, t, 0.00370026744188713119232403)
    u = mla(u, t, -0.00889896195887655491740809)
    u = mla(u, t, 0.016599329773529201970117)
    u = mla(u, t, -0.0254517624932312641616861)
    u = mla(u, t, 0.0337852580001353069993897)
    u = mla(u, t, -0.0407629191276836500001934)
    u = mla(u, t, 0.0466667150077840625632675)
    u = mla(u, t, -0.0523674852303482457616113)
    u = mla(u, t, 0.0587666392926673580854313)
    u = mla(u, t, -0.0666573579361080525984562)
    u = mla(u, t, 0.0769219538311769618355029)
    u = mla(u, t, -0.090908995008245008229153)
    u = mla(u, t, 0.111111105648261418443745)
    u = mla(u, t, -0.14285714266771329383765)
    u = mla(u, t, 0.199999999996591265594148)
    u = mla(u, t, -0.333333333333311110369124)

    t = s + s*(t*u)

    (q & 1) != 0 && (t = 1.570796326794896557998982 - t)
    (q & 2) != 0 && (t = -t)

    return t
end

function xatan2_u1(y::Float64, x::Float64)
    d = atan2k_u1(Double2(abs(y), 0.0), Double2(x, 0.0))
    r = d.x + d.y

    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI/2)) : 0)
    end
    if isinf(y)
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI*1/4)) : 0)
    end
    if y == 0
        r = sign(x) == -1 ? M_PI : 0
    end
    return isnan(x) || isnan(y) ? NaN : flipsign(r, y)
end

function xasin_u1(d::Float64)
    d2 = atan2k_u1(Double2(abs(d), 0.0), ddsqrt_d2_d2(ddmul_d2_d2_d2(ddadd_d2_d_d(1.0, d), ddadd_d2_d_d(1.0,-d))))
    r = d2.x + d2.y
    abs(d) == 1 && (r = 1.570796326794896557998982)
    return flipsign(r, d)
end

function xacos_u1(d::Float64)
    d2 = atan2k_u1(ddsqrt_d2_d2(ddmul_d2_d2_d2(ddadd_d2_d_d(1.0, d), ddadd_d2_d_d(1.0,-d))), Double2(abs(d), 0.0))
    d2 = ddscale_d2_d2_d(d2, sign(d))
    abs(d) == 1 && (d2 = Double2(0.0, 0.0))
    d < 0 && (d2 = ddadd_d2_d2_d2(Double2(3.141592653589793116, 1.2246467991473532072e-16), d2))
    return d2.x + d2.y
end

function xatan_u1(d::Float64)
    d2 = atan2k_u1(Double2(abs(d), 0.0), Double2(1.0, 0.0))
    r = d2.x + d2.y
    isinf(d) && (r = 1.570796326794896557998982)
    return flipsign(r, d)
end

function xsin(d::Float64)
    q = xrint(d*M_1_PI)

    d = mla(q, -PI4_A*4, d)
    d = mla(q, -PI4_B*4, d)
    d = mla(q, -PI4_C*4, d)
    d = mla(q, -PI4_D*4, d)

    s = d*d

    (q & 1) != 0 && (d = -d)

    u = -7.97255955009037868891952e-18
    u = mla(u, s, 2.81009972710863200091251e-15)
    u = mla(u, s, -7.64712219118158833288484e-13)
    u = mla(u, s, 1.60590430605664501629054e-10)
    u = mla(u, s, -2.50521083763502045810755e-08)
    u = mla(u, s, 2.75573192239198747630416e-06)
    u = mla(u, s, -0.000198412698412696162806809)
    u = mla(u, s, 0.00833333333333332974823815)
    u = mla(u, s, -0.166666666666666657414808)

    u = mla(s, u*d, d)
    return u
end

function xsin_u1(d::Float64)
    q = xrint(d*M_1_PI)

    s = ddadd2_d2_d_d(d, q * (-PI4_A*4))
    s = ddadd2_d2_d2_d(s, q * (-PI4_B*4))
    s = ddadd2_d2_d2_d(s, q * (-PI4_C*4))
    s = ddadd2_d2_d2_d(s, q * (-PI4_D*4))

    t = s
    s = ddsqu_d2_d2(s)

    u = 2.72052416138529567917983e-15
    u = mla(u, s.x, -7.6429259411395447190023e-13)
    u = mla(u, s.x, 1.60589370117277896211623e-10)
    u = mla(u, s.x, -2.5052106814843123359368e-08)
    u = mla(u, s.x, 2.75573192104428224777379e-06)
    u = mla(u, s.x, -0.000198412698412046454654947)
    u = mla(u, s.x, 0.00833333333333318056201922)

    x = ddadd_d2_d_d2(1.0, ddmul_d2_d2_d2(ddadd_d2_d_d(-0.166666666666666657414808, u * s.x), s))

    x = ddmul_d2_d2_d2(t, x)
    u = x.x + x.y

    (q & 1) != 0 && (u = -u)

    return u
end

function xcos(d::Float64)
    q = 1 + 2*xrint(d * M_1_PI - 0.5)

    d = mla(q, -PI4_A*2, d)
    d = mla(q, -PI4_B*2, d)
    d = mla(q, -PI4_C*2, d)
    d = mla(q, -PI4_D*2, d)

    s = d*d

    (q & 2) == 0 && (d = -d)

    u = -7.97255955009037868891952e-18
    u = mla(u, s, 2.81009972710863200091251e-15)
    u = mla(u, s, -7.64712219118158833288484e-13)
    u = mla(u, s, 1.60590430605664501629054e-10)
    u = mla(u, s, -2.50521083763502045810755e-08)
    u = mla(u, s, 2.75573192239198747630416e-06)
    u = mla(u, s, -0.000198412698412696162806809)
    u = mla(u, s, 0.00833333333333332974823815)
    u = mla(u, s, -0.166666666666666657414808)

    u = mla(s, u * d, d)

    return u
end

function xcos_u1(d::Float64)
    d = abs(d)

    q = mla(2, xrint(d * M_1_PI - 0.5), 1)

    s = ddadd2_d2_d_d(d, q * (-PI4_A*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_B*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_C*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_D*2))

    t = s
    s = ddsqu_d2_d2(s)

    u = 2.72052416138529567917983e-15
    u = mla(u, s.x, -7.6429259411395447190023e-13)
    u = mla(u, s.x, 1.60589370117277896211623e-10)
    u = mla(u, s.x, -2.5052106814843123359368e-08)
    u = mla(u, s.x, 2.75573192104428224777379e-06)
    u = mla(u, s.x, -0.000198412698412046454654947)
    u = mla(u, s.x, 0.00833333333333318056201922)

    x = ddadd_d2_d_d2(1.0, ddmul_d2_d2_d2(ddadd_d2_d_d(-0.166666666666666657414808, u * s.x), s))
    x = ddmul_d2_d2_d2(t, x)

    u = x.x + x.y
    (q & 2) == 0 && (u = -u)

    return u
end

function xsincos(d::Float64)
    q = xrint(d*(2*M_1_PI))
    s = d

    s = mla(-q, PI4_A*2, s)
    s = mla(-q, PI4_B*2, s)
    s = mla(-q, PI4_C*2, s)
    s = mla(-q, PI4_D*2, s)

    t = s
    s = s*s

    u = 1.58938307283228937328511e-10
    u = mla(u, s, -2.50506943502539773349318e-08)
    u = mla(u, s, 2.75573131776846360512547e-06)
    u = mla(u, s, -0.000198412698278911770864914)
    u = mla(u, s, 0.0083333333333191845961746)
    u = mla(u, s, -0.166666666666666130709393)
    u = u * s * t

    rx = t + u

    u = -1.13615350239097429531523e-11
    u = mla(u, s, 2.08757471207040055479366e-09)
    u = mla(u, s, -2.75573144028847567498567e-07)
    u = mla(u, s, 2.48015872890001867311915e-05)
    u = mla(u, s, -0.00138888888888714019282329)
    u = mla(u, s, 0.0416666666666665519592062)
    u = mla(u, s, -0.5)

    ry = u * s + 1

    (q & 1) != 0 && (s = ry; ry = rx; rx = s)
    (q & 2) != 0 && (rx = -rx)
    ((q+1) & 2) != 0 && (ry = -ry)

    isinf(d) && (rx = ry = NaN)

    return Double2(rx,ry)
end

function xsincos_u1(d::Float64)
    q = xrint(d*(2*M_1_PI))

    s = ddadd2_d2_d_d(d, q * (-PI4_A*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_B*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_C*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_D*2))

    t = s
    s = ddsqu_d2_d2(s)
    sx = s.x + s.y

    u = 1.58938307283228937328511e-10
    u = mla(u, sx, -2.50506943502539773349318e-08)
    u = mla(u, sx, 2.75573131776846360512547e-06)
    u = mla(u, sx, -0.000198412698278911770864914)
    u = mla(u, sx, 0.0083333333333191845961746)
    u = mla(u, sx, -0.166666666666666130709393)

    u *= sx * t.x

    x = ddadd_d2_d2_d(t, u)
    rx = x.x + x.y

    u = -1.13615350239097429531523e-11
    u = mla(u, sx, 2.08757471207040055479366e-09)
    u = mla(u, sx, -2.75573144028847567498567e-07)
    u = mla(u, sx, 2.48015872890001867311915e-05)
    u = mla(u, sx, -0.00138888888888714019282329)
    u = mla(u, sx, 0.0416666666666665519592062)
    u = mla(u, sx, -0.5)

    x = ddadd_d2_d_d2(1.0, ddmul_d2_d_d(sx, u))
    ry = x.x + x.y

    (q & 1) != 0 && (u = ry; ry = rx; rx = u)
    (q & 2) != 0 && (rx = -rx)
    ((q+1) & 2) != 0 && (ry = -ry)

    isinf(d) && (rx = ry = NaN)

    return Double2(rx,ry)
end

function xtan(d::Float64)
    q = xrint(d * (2 * M_1_PI))

    x = mla(q, -PI4_A*2, d)
    x = mla(q, -PI4_B*2, x)
    x = mla(q, -PI4_C*2, x)
    x = mla(q, -PI4_D*2, x)

    s = x*x

    (q & 1) != 0 && (x = -x)

    u = 1.01419718511083373224408e-05;
    u = mla(u, s, -2.59519791585924697698614e-05)
    u = mla(u, s, 5.23388081915899855325186e-05)
    u = mla(u, s, -3.05033014433946488225616e-05)
    u = mla(u, s, 7.14707504084242744267497e-05)
    u = mla(u, s, 8.09674518280159187045078e-05)
    u = mla(u, s, 0.000244884931879331847054404)
    u = mla(u, s, 0.000588505168743587154904506)
    u = mla(u, s, 0.00145612788922812427978848)
    u = mla(u, s, 0.00359208743836906619142924)
    u = mla(u, s, 0.00886323944362401618113356)
    u = mla(u, s, 0.0218694882853846389592078)
    u = mla(u, s, 0.0539682539781298417636002)
    u = mla(u, s, 0.133333333333125941821962)
    u = mla(u, s, 0.333333333333334980164153)

    u = mla(s, u * x, x)

    (q & 1) != 0 && (u = 1.0/u)
    isinf(d) && (u = NaN)

    return u
end

function xtan_u1(d::Float64)
    q = xrint(d * M_2_PI)

    s = ddadd2_d2_d_d(d, q * (-PI4_A*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_B*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_C*2))
    s = ddadd2_d2_d2_d(s, q * (-PI4_D*2))

    (q & 1) != 0 && (s = ddneg_d2_d2(s))

    t = s
    s = ddsqu_d2_d2(s)

    u = 1.01419718511083373224408e-05
    u = mla(u, s.x, -2.59519791585924697698614e-05)
    u = mla(u, s.x, 5.23388081915899855325186e-05)
    u = mla(u, s.x, -3.05033014433946488225616e-05)
    u = mla(u, s.x, 7.14707504084242744267497e-05)
    u = mla(u, s.x, 8.09674518280159187045078e-05)
    u = mla(u, s.x, 0.000244884931879331847054404)
    u = mla(u, s.x, 0.000588505168743587154904506)
    u = mla(u, s.x, 0.00145612788922812427978848)
    u = mla(u, s.x, 0.00359208743836906619142924)
    u = mla(u, s.x, 0.00886323944362401618113356)
    u = mla(u, s.x, 0.0218694882853846389592078)
    u = mla(u, s.x, 0.0539682539781298417636002)
    u = mla(u, s.x, 0.133333333333125941821962)

    x = ddadd_d2_d_d2(1.0, ddmul_d2_d2_d2(ddadd_d2_d_d(0.333333333333334980164153, u * s.x), s))
    x = ddmul_d2_d2_d2(t, x)

    (q & 1) != 0 && (x = ddrec_d2_d2(x))

    u = x.x + x.y

    return u
end

function xlog(d::Float64)
    e = ilogbp1(d*0.7071)
    m = ldexpk(d,-e)

    x = (m-1)/(m+1)
    x2 = x*x

    t = 0.148197055177935105296783;
    t = mla(t, x2, 0.153108178020442575739679)
    t = mla(t, x2, 0.181837339521549679055568)
    t = mla(t, x2, 0.22222194152736701733275)
    t = mla(t, x2, 0.285714288030134544449368)
    t = mla(t, x2, 0.399999999989941956712869)
    t = mla(t, x2, 0.666666666666685503450651)
    t = mla(t, x2, 2)

    x = x * t + 0.693147180559945286226764 * e

    isinf(d) && (x = Inf)
    d < 0 && (x = NaN)
    d == 0 && (x = -Inf)

    return x
end

function xexp(d::Float64)
    q = xrint(d*R_LN2)

    s = mla(q,-L2U, d)
    s = mla(q,-L2L, s)

    u = 2.08860621107283687536341e-09;
    u = mla(u, s, 2.51112930892876518610661e-08)
    u = mla(u, s, 2.75573911234900471893338e-07)
    u = mla(u, s, 2.75572362911928827629423e-06)
    u = mla(u, s, 2.4801587159235472998791e-05)
    u = mla(u, s, 0.000198412698960509205564975)
    u = mla(u, s, 0.00138888888889774492207962)
    u = mla(u, s, 0.00833333333331652721664984)
    u = mla(u, s, 0.0416666666666665047591422)
    u = mla(u, s, 0.166666666666666851703837)
    u = mla(u, s, 0.5)

    u = s * s * u + s + 1
    u = ldexpk(u, q)

    d == -Inf && (u = 0.0)
    
    return u
end

function xlog_u1(d::Float64)
    s = logk(d)
    x = s.x + s.y

    isinf(d) && (x = Inf)
    d < 0    && (x = NaN)
    d == 0   && (x = -Inf)
    return x
end

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
    y = abs(x);
    d = logk2(ddadd_d2_d2_d(ddsqrt_d2_d2(ddadd2_d2_d2_d(ddmul_d2_d_d(y, y),  1.0)), y))
    y = d.x + d.y

    y = isinf(x) || isnan(y) ? Inf : y
    y = flipsign(y, x)
    y = isnan(x) ? NaN : y

    return y
end

function xacosh(x::Float64)
    d = logk2(ddadd2_d2_d2_d(ddsqrt_d2_d2(ddadd2_d2_d2_d(ddmul_d2_d_d(x, x), -1)), x))
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

function xcbrt(d::Float64) # max error : 2 ulps
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

function xexp2(a::Float64)
    u = expk(ddmul_d2_d2_d(Double2(0.69314718055994528623, 2.3190468138462995584e-17), a))
    a > 1023   && (u = Inf)
    a == -Inf && (u = 0.0)
  return u
end

function xexp10(a::Float64)
    u = expk(ddmul_d2_d2_d(Double2(2.3025850929940459011, -2.1707562233822493508e-16), a))
    a > 308   && (u = Inf)
    a == -Inf && (u = 0.0)
    return u
end

function xexpm1(a::Float64)
    d = ddadd2_d2_d2_d(expk2(Double2(a, 0.0)), -1.0)
    x = d.x + d.y
    a > 700 && (x = Inf)
    a < -0.36043653389117156089696070315825181539851971360337e+2 && (x = -1.0)
    return x
end

function xlog10(a::Float64)
    d = ddmul_d2_d2_d2(logk(a), Double2(0.43429448190325176116, 6.6494347733425473126e-17))
    x = d.x + d.y
    isinf(a) && (x = Inf)
    a < 0    && (x = NaN)
    a == 0   && (x = -Inf)
    return x
end

function xlog1p(a::Float64)
    d = logk2(ddadd2_d2_d_d(a, 1.0))
    x = d.x + d.y
    isinf(a) && (x = Inf)
    a < -1   && (x = NaN)
    a == -1  && (x = -Inf)
    return x
end
end