
# utility functions

# remove and just use reinterpret instead of the two function below (todo)
double_to_raw_long_bits(d::Float64) = reinterpret(Int64, d)
long_bits_to_double(i::Int64) = reinterpret(Float64, i)

# this seems to emit better native code than Base.sign
sign(d::Float64) =  flipsign(1.0, d)

mla(x::Number, y::Number, z::Number) = muladd(x, y, z)

xrint(x::Float64) = unsafe_trunc(Int32, round(x)) # in sleef, but this is a limited way to truncate since Int32 (fix)

pow2i(q::Int32) = long_bits_to_double(((q + Int32(0x3ff)) % Int64) << 52)


# private math functions

function ldexpk(x::Float64, q::Int32)
    m = q >> 31
    m = (((m + q) >> 9) - m) << 7
    q = q - (m << 2)
    m += Int32(0x3ff)
    m = m < 0 ? Int32(0) : m
    m = m > 0x7ff ? Int32(0x7ff) : m
    u = long_bits_to_double((m % Int64) << 52)
    x = x * u * u * u * u
    u = long_bits_to_double(((q + Int32(0x3ff)) % Int64) << 52)
    return x * u
end

function ilogbp1(d::Float64)
    m = d < 4.9090934652977266e-91
    d = m ? 2.037035976334486e90 * d : d
    q = ((double_to_raw_long_bits(d) >> 52) & 0x7ff) 
    q = m ? q - (300 + 0x03fe) : q - 0x03fe
    return q % Int32
end

function atan2k(y::Float64, x::Float64)
    q = 0
    if x < 0
        x = -x
        q = -2
    end
    if y > x
        t = x
        x = y
        y = -t
        q += 1
    end
    s = y/x
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

    t = u*t*s + s
    t = q*(M_PI/2) + t
    return t
end

function atan2k_u1(y::Double2, x::Double2)
    q = 0
    if x.x < 0
        x = Double2(-x.x,-x.y)
        q = -2
    end
    if y.x > x.x
        t = x
        x = y
        y = Double2(-t.x,-t.y)
        q += 1
    end
    s = dddiv_d2_d2_d2(y, x)
    t = ddsqu_d2_d2(s)
    t = ddnormalize_d2_d2(t)

    u = 1.06298484191448746607415e-05
    u = mla(u, t.x, -0.000125620649967286867384336)
    u = mla(u, t.x, 0.00070557664296393412389774)
    u = mla(u, t.x, -0.00251865614498713360352999)
    u = mla(u, t.x, 0.00646262899036991172313504)
    u = mla(u, t.x, -0.0128281333663399031014274)
    u = mla(u, t.x, 0.0208024799924145797902497)
    u = mla(u, t.x, -0.0289002344784740315686289)
    u = mla(u, t.x, 0.0359785005035104590853656)
    u = mla(u, t.x, -0.041848579703592507506027)
    u = mla(u, t.x, 0.0470843011653283988193763)
    u = mla(u, t.x, -0.0524914210588448421068719)
    u = mla(u, t.x, 0.0587946590969581003860434)
    u = mla(u, t.x, -0.0666620884778795497194182)
    u = mla(u, t.x, 0.0769225330296203768654095)
    u = mla(u, t.x, -0.0909090442773387574781907)
    u = mla(u, t.x, 0.111111108376896236538123)
    u = mla(u, t.x, -0.142857142756268568062339)
    u = mla(u, t.x, 0.199999999997977351284817)
    u = mla(u, t.x, -0.333333333333317605173818)

    t = ddmul_d2_d2_d(t, u)
    t = ddmul_d2_d2_d2(s, ddadd_d2_d_d2(1.0, t))
    t = ddadd2_d2_d2_d2(ddmul_d2_d2_d(Double2(1.570796326794896557998982, 6.12323399573676603586882e-17), Float64(q)), t)
    return t
end

function logk(d::Float64)
    e = ilogbp1(d*0.7071)
    m = ldexpk(d,-e)

    x = dddiv_d2_d2_d2(ddadd2_d2_d_d(-1.0, m), ddadd2_d2_d_d(1.0, m))
    x2 = ddsqu_d2_d2(x)

    t = 0.134601987501262130076155
    t = mla(t, x2.x, 0.132248509032032670243288)
    t = mla(t, x2.x, 0.153883458318096079652524)
    t = mla(t, x2.x, 0.181817427573705403298686)
    t = mla(t, x2.x, 0.222222231326187414840781)
    t = mla(t, x2.x, 0.285714285651261412873718)
    t = mla(t, x2.x, 0.400000000000222439910458)
    t = mla(t, x2.x, 0.666666666666666371239645)

    return ddadd2_d2_d2_d2(ddmul_d2_d2_d(Double2(0.693147180559945286226764, 2.319046813846299558417771e-17), Float64(e)),
            ddadd2_d2_d2_d2(ddscale_d2_d2_d(x, 2.0), ddmul_d2_d2_d(ddmul_d2_d2_d2(x2, x), t)))
end

function expk(d::Double2)
    q = xrint((d.x + d.y) * R_LN2)

    s = ddadd2_d2_d2_d(d, q*-L2U)
    s = ddadd2_d2_d2_d(s, q*-L2L)

    s = ddnormalize_d2_d2(s)

    u = 2.51069683420950419527139e-08
    u = mla(u, s.x, 2.76286166770270649116855e-07)
    u = mla(u, s.x, 2.75572496725023574143864e-06)
    u = mla(u, s.x, 2.48014973989819794114153e-05)
    u = mla(u, s.x, 0.000198412698809069797676111)
    u = mla(u, s.x, 0.0013888888939977128960529)
    u = mla(u, s.x, 0.00833333333332371417601081)
    u = mla(u, s.x, 0.0416666666665409524128449)
    u = mla(u, s.x, 0.166666666666666740681535)
    u = mla(u, s.x, 0.500000000000000999200722)

    t = ddadd_d2_d2_d2(s, ddmul_d2_d2_d(ddsqu_d2_d2(s), u))

    t = ddadd_d2_d_d2(1.0, t)
    return ldexpk(t.x + t.y, q)
end

function expk2(d::Double2)
    q = xrint((d.x + d.y)*R_LN2)

    s = ddadd2_d2_d2_d(d, q*-L2U)
    s = ddadd2_d2_d2_d(s, q*-L2L)

    u = 2.51069683420950419527139e-08
    u = mla(u, s.x, 2.76286166770270649116855e-07)
    u = mla(u, s.x, 2.75572496725023574143864e-06)
    u = mla(u, s.x, 2.48014973989819794114153e-05)
    u = mla(u, s.x, 0.000198412698809069797676111)
    u = mla(u, s.x, 0.0013888888939977128960529)
    u = mla(u, s.x, 0.00833333333332371417601081)
    u = mla(u, s.x, 0.0416666666665409524128449)
    u = mla(u, s.x, 0.166666666666666740681535)
    u = mla(u, s.x, 0.500000000000000999200722)

    t = ddadd_d2_d2_d2(s, ddmul_d2_d2_d(ddsqu_d2_d2(s), u))

    t = ddadd_d2_d_d2(1.0, t)
    return ddscale_d2_d2_d(t, pow2i(q))
end

function logk2(d::Double2)
    e = ilogbp1(d.x * 0.7071)
    m = ddscale_d2_d2_d(d, pow2i(-e))

    x = dddiv_d2_d2_d2(ddadd2_d2_d2_d(m, -1.0), ddadd2_d2_d2_d(m, 1.0))
    x2 = ddsqu_d2_d2(x)

    t = 0.134601987501262130076155
    t = mla(t, x2.x, 0.132248509032032670243288)
    t = mla(t, x2.x, 0.153883458318096079652524)
    t = mla(t, x2.x, 0.181817427573705403298686)
    t = mla(t, x2.x, 0.222222231326187414840781)
    t = mla(t, x2.x, 0.285714285651261412873718)
    t = mla(t, x2.x, 0.400000000000222439910458)
    t = mla(t, x2.x, 0.666666666666666371239645)

    return ddadd2_d2_d2_d2(ddmul_d2_d2_d(Double2(0.693147180559945286226764, 2.319046813846299558417771e-17), Float64(e)),
            ddadd2_d2_d2_d2(ddscale_d2_d2_d(x, 2.0), ddmul_d2_d2_d(ddmul_d2_d2_d2(x2, x), t)))
end