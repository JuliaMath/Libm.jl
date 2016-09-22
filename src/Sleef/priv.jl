# private math functions

"""
A helper function for `ldexpk`

First note that `r = (q >> n) << n` clears the lowest n bits of q, i.e. returns 2^n where n is the
largest integer such that q >= 2^n

For numbers q less than 2^m the following code does the same as the above snippet
    `r = ( (q>>v + q) >> n - q>>v ) << n`
For numbers larger than or equal to 2^v this subtracts 2^n from q for q>>n times.

The function returns q(input) := q(output) + offset*r

In the code for ldexpk we actually use
    `m = ( (m>>n + m) >> n -  m>>m ) << (n-2)`.
So that x has to be multplied by u four times `x = x*u*u*u*u` to put the value  of the offset
exponent amount back in.
"""
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    return m, q
end
@inline split_exponent(::Type{Float64}, q::Int) = _split_exponent(q, UInt(9), UInt(31), UInt(2))
@inline split_exponent(::Type{Float32}, q::Int) = _split_exponent(q, UInt(6), UInt(31), UInt(2))

"""
    ldexpk(x::FloatTypes, n::Int) -> FloatTypes

Computes `x \times 2^n`
"""
@inline function ldexpk{T<:FloatTypes}(x::T, q::Int)
    bias = exponent_bias(T)
    emax = exponent_max(T)
    m, q = split_exponent(T,q)
    m += bias
    m = m < 0 ? 0 : m
    m = m > emax ? emax : m
    q += bias
    u = integer2float(T, m)
    x = x*u*u*u*u
    u = integer2float(T, q)
    return x*u
end

"""
    exponent = ilogbp1(x)

Returns the integral part of the logarithm of `|x|`, using 2 as base for the logarithm; in other
words this returns the binary exponent of `x` so that
    x = significand \times 2^exponenet
where `significand \in [0.5, 1)`
"""

# The following define threshold values for `ilogbp1`
real_cut_offset(::Type{Float64}) = 300
real_cut_offset(::Type{Float32}) = 64

# 2^-real_cut_offset
real_cut_min(::Type{Float64}) = 4.9090934652977266e-91
real_cut_min(::Type{Float32}) = 5.421010862427522f-20

# 2^real_cut_offset
real_cut_max(::Type{Float64}) = 2.037035976334486e90
real_cut_max(::Type{Float32}) = 1.8446744073709552f19

@inline function ilogbp1{T<:FloatTypes}(d::T)
    m = d < real_cut_min(T)
    d = m ? real_cut_max(T) * d : d
    q = float2integer(d) & exponent_max(T)
    q = m ? q - (real_cut_offset(T) + exponent_bias(T) - 1) : q - (exponent_bias(T) - 1) # we subtract 1 since we want 2^q
    return q
end


let
global atan2k
const c19d = -1.88796008463073496563746e-05
const c18d =  0.000209850076645816976906797
const c17d = -0.00110611831486672482563471
const c16d =  0.00370026744188713119232403
const c15d = -0.00889896195887655491740809
const c14d =  0.016599329773529201970117
const c13d = -0.0254517624932312641616861
const c12d =  0.0337852580001353069993897
const c11d = -0.0407629191276836500001934
const c10d =  0.0466667150077840625632675
const c9d  = -0.0523674852303482457616113
const c8d  =  0.0587666392926673580854313
const c7d  = -0.0666573579361080525984562
const c6d  =  0.0769219538311769618355029
const c5d  = -0.090908995008245008229153
const c4d  =  0.111111105648261418443745
const c3d  = -0.14285714266771329383765
const c2d  =  0.199999999996591265594148
const c1d  = -0.333333333333311110369124

const c8f =  0.00282363896258175373077393f0
const c7f = -0.0159569028764963150024414f0
const c6f =  0.0425049886107444763183594f0
const c5f = -0.0748900920152664184570312f0
const c4f =  0.106347933411598205566406f0
const c3f = -0.142027363181114196777344f0
const c2f =  0.199926957488059997558594f0
const c1f = -0.333331018686294555664062f0

global @inline _atan2k(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d
global @inline _atan2k(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f

@inline function atan2k{T<:FloatTypes}(y::T, x::T)
    q = 0
    if x < 0
        x = -x
        q = -2
    end
    if y > x
        t = x; x = y;
        y = -t
        q += 1
    end
    s = y/x
    t = s*s
    u = _atan2k(t)
    t = u*t*s + s
    return q*T(MPI/2) + t
end
end

let
global atan2k_u1
const c20 =  1.06298484191448746607415e-05
const c19 = -0.000125620649967286867384336
const c18 =  0.00070557664296393412389774
const c17 = -0.00251865614498713360352999
const c16 =  0.00646262899036991172313504
const c15 = -0.0128281333663399031014274
const c14 =  0.0208024799924145797902497
const c13 = -0.0289002344784740315686289
const c12 =  0.0359785005035104590853656
const c11 = -0.041848579703592507506027
const c10 =  0.0470843011653283988193763
const c9  = -0.0524914210588448421068719
const c8  =  0.0587946590969581003860434
const c7  = -0.0666620884778795497194182
const c6  =  0.0769225330296203768654095
const c5  = -0.0909090442773387574781907
const c4  =  0.111111108376896236538123
const c3  = -0.142857142756268568062339
const c2  =  0.199999999997977351284817
const c1  = -0.333333333333317605173818

@inline function atan2k_u1{T<:Float64}(y::Double{T}, x::Double{T})
    q = 0
    if x.hi < 0
        x = Double(-x.hi,-x.lo)
        q = -2
    end
    if y.hi > x.hi
        t = x; x = y
        y = Double(-t.hi,-t.lo)
        q += 1
    end
    s = dddiv(y, x)
    t = ddsqu(s)
    t = ddnormalize(t)
    u = @horner t.hi c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20
    t = ddmul(t, u)
    t = ddmul(s, ddadd(T(1.0), t))
    return ddadd2(ddmul(Double(1.570796326794896557998982, 6.12323399573676603586882e-17), T(q)), t)
end
end

let
global logk
const c8 = 0.134601987501262130076155
const c7 = 0.132248509032032670243288
const c6 = 0.153883458318096079652524
const c5 = 0.181817427573705403298686
const c4 = 0.222222231326187414840781
const c3 = 0.285714285651261412873718
const c2 = 0.400000000000222439910458
const c1 = 0.666666666666666371239645

@inline function logk{T<:Float64}(d::T)
    e = ilogbp1(d*0.7071)
    m = ldexpk(d,-e)
    x = dddiv(ddadd2(-1.0, m), ddadd2(1.0, m))
    x2 = ddsqu(x)
    t = @horner x2.hi c1 c2 c3 c4 c5 c6 c7 c8
    return ddadd2(ddmul(Double(0.693147180559945286226764, 2.319046813846299558417771e-17), T(e)),
            ddadd2(ddscale(x, T(2.0)), ddmul(ddmul(x2, x), t)))
end
end

let
global expk
const c10 = 2.51069683420950419527139e-08
const c9 = 2.76286166770270649116855e-07
const c8 = 2.75572496725023574143864e-06
const c7 = 2.48014973989819794114153e-05
const c6 = 0.000198412698809069797676111
const c5 = 0.0013888888939977128960529
const c4 = 0.00833333333332371417601081
const c3 = 0.0416666666665409524128449
const c2 = 0.166666666666666740681535
const c1 = 0.500000000000000999200722

@inline function expk{T<:Float64}(d::Double{T})
    q = xrint((d.hi + d.lo)*LOG2E)
    s = ddadd2(d, q*-LN2U(T))
    s = ddadd2(s, q*-LN2L(T))
    s = ddnormalize(s)
    u = @horner s.hi c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
    t = ddadd(s, ddmul(ddsqu(s), u))
    t = ddadd(one(T), t)
    return ldexpk(t.hi + t.lo, q)
end
end

let
global expk2
const c10 = 2.51069683420950419527139e-08
const c9 = 2.76286166770270649116855e-07
const c8 = 2.75572496725023574143864e-06
const c7 = 2.48014973989819794114153e-05
const c6 = 0.000198412698809069797676111
const c5 = 0.0013888888939977128960529
const c4 = 0.00833333333332371417601081
const c3 = 0.0416666666665409524128449
const c2 = 0.166666666666666740681535
const c1 = 0.500000000000000999200722

@inline function expk2{T<:Float64}(d::Double{T})
    q = xrint((d.hi + d.lo)*T(LOG2E))
    s = ddadd2(d, q * -LN2U(T))
    s = ddadd2(s, q * -LN2L(T))
    u = @horner s.hi c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
    t = ddadd(s, ddmul(ddsqu(s), u))
    t = ddadd(T(1.0), t)
    return ddscale(t, pow2i(T, q))
end
end

let
global logk2
const c8 = 0.134601987501262130076155
const c7 = 0.132248509032032670243288
const c6 = 0.153883458318096079652524
const c5 = 0.181817427573705403298686
const c4 = 0.222222231326187414840781
const c3 = 0.285714285651261412873718
const c2 = 0.400000000000222439910458
const c1 = 0.666666666666666371239645

@inline function logk2{T<:Float64}(d::Double{T})
    e = ilogbp1(d.hi * 0.7071)
    m = ddscale(d, pow2i(T, -e))
    x = dddiv(ddadd2(m, T(-1.0)), ddadd2(m, T(1.0)))
    x2 = ddsqu(x)
    t = @horner x2.hi c1 c2 c3 c4 c5 c6 c7 c8
    return ddadd2(ddmul(Double(0.693147180559945286226764, 2.319046813846299558417771e-17), T(e)),
            ddadd2(ddscale(x, T(2.0)), ddmul(ddmul(x2, x), t)))
end
end