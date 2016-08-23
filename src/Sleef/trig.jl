
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

function xatan2(y::Float64, x::Float64)
    r = atan2k(abs(y), x)

    r = flipsign(r, x)
    if isinf(x) || x == 0
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI/2)) : 0.0)
    end
    if isinf(y)
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI/4)) : 0.0)
    end
    if y == 0
        r = (sign(x) == -1 ? M_PI/1 : 0.0)
    end
    return isnan(x) || isnan(y) ? NaN : flipsign(r, y)
end

xasin(d::Float64) = flipsign(atan2k(abs(d), sqrt((1+d)*(1-d))), d)

xacos(d::Float64) = flipsign(atan2k(sqrt((1+d)*(1-d)), abs(d)), d) + (d < 0 ? M_PI/1 : 0.0)

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
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI/2)) : 0.0)
    end
    if isinf(y)
        r = M_PI/2 - (isinf(x) ? (sign(x) * (M_PI/4)) : 0.0)
    end
    if y == 0
        r = sign(x) == -1 ? M_PI/1 : 0.0
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
