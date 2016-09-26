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

# The following define threshold values for `ilogbp1`
real_cut_offset(::Type{Float64}) = 300
real_cut_offset(::Type{Float32}) = 64

# 2^-real_cut_offset
real_cut_min(::Type{Float64}) = 4.9090934652977266e-91
real_cut_min(::Type{Float32}) = 5.421010862427522f-20

# 2^real_cut_offset
real_cut_max(::Type{Float64}) = 2.037035976334486e90
real_cut_max(::Type{Float32}) = 1.8446744073709552f19

"""
    ilogbp1(x::FloatTypes) -> Int

Returns the integral part of the logarithm of `|x|`, using 2 as base for the logarithm; in other
words this returns the binary exponent of `x` so that
    x = significand \times 2^exponenet
where `significand \in [0.5, 1)`
"""
@inline function ilogbp1{T<:FloatTypes}(d::T)
    m = d < real_cut_min(T)
    d = m ? real_cut_max(T) * d : d
    q = float2integer(d) & exponent_max(T)
    q = m ? q - (real_cut_offset(T) + exponent_bias(T) - 1) : q - (exponent_bias(T) - 1) # we subtract 1 since we want 2^q
    return q
end


let
const c20d =  1.06298484191448746607415e-05
const c19d = -0.000125620649967286867384336
const c18d =  0.00070557664296393412389774
const c17d = -0.00251865614498713360352999
const c16d =  0.00646262899036991172313504
const c15d = -0.0128281333663399031014274
const c14d =  0.0208024799924145797902497
const c13d = -0.0289002344784740315686289
const c12d =  0.0359785005035104590853656
const c11d = -0.041848579703592507506027
const c10d =  0.0470843011653283988193763
const c9d  = -0.0524914210588448421068719
const c8d  =  0.0587946590969581003860434
const c7d  = -0.0666620884778795497194182
const c6d  =  0.0769225330296203768654095
const c5d  = -0.0909090442773387574781907
const c4d  =  0.111111108376896236538123
const c3d  = -0.142857142756268568062339
const c2d  =  0.199999999997977351284817
const c1d =  -0.333333333333317605173818

const c9f = -0.00176397908944636583328247f0
const c8f =  0.0107900900766253471374512f0
const c7f = -0.0309564601629972457885742f0
const c6f =  0.0577365085482597351074219f0
const c5f = -0.0838950723409652709960938f0
const c4f =  0.109463557600975036621094f0
const c3f = -0.142626821994781494140625f0
const c2f =  0.199983194470405578613281f0
const c1f = -0.333332866430282592773438f0

global @inline _atan2k(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline _atan2k(x::Float32) = @horner x c1f c2f c3f c4f c5f c6f c7f c8f c9f

global @inline function atan2k{T<:FloatTypes}(y::T, x::T)
    q = 0
    if x < 0
        x = -x
        q = -2
    end
    if y > x
        t = x; x = y
        y = -t
        q += 1
    end
    s = y/x
    t = s*s
    u =_atan2k(t)
    t = u*t*s + s
    return q*T(MPI2) + t
end

global @inline _atan2k_u1(x::Double{Float64}) = @horner x.hi c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d c16d c17d c18d c19d c20d
global @inline _atan2k_u1(x::Double{Float32}) = ddadd(c1f, x.hi*(@horner x.hi c2f c3f c4f c5f c6f c7f c8f c9f))

global @inline function atan2k_u1{T<:FloatTypes}(y::Double{T}, x::Double{T})
    q = 0
    if x.hi < 0
        x = -x
        q = -2
    end
    if y.hi > x.hi
        t = x; x = y
        y = -t
        q += 1
    end
    s = dddiv(y, x)
    t = ddsqu(s)
    t = normalize(t)
    u =_atan2k_u1(t)
    t = ddmul(s, ddadd(T(1), ddmul(t, u)))
    return ddadd2(ddmul(T(q), MDPI2(T)), t)
end
end

let
const c10d = 2.51069683420950419527139e-08
const c9d  = 2.76286166770270649116855e-07
const c8d  = 2.75572496725023574143864e-06
const c7d  = 2.48014973989819794114153e-05
const c6d  = 0.000198412698809069797676111
const c5d  = 0.001388888893997712896052900
const c4d  = 0.008333333333323714176010810
const c3d  = 0.041666666666540952412844900
const c2d  = 0.166666666666666740681535000
const c1d  = 0.500000000000000999200722000

const c5f = 0.00136324646882712841033936f0
const c4f = 0.00836596917361021041870117f0
const c3f = 0.04167108237743377685546880f0
const c2f = 0.16666552424430847167968800f0
const c1f = 0.49999985098838806152343800f0

global @inline _expk(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d
global @inline _expk(x::Float32) = @horner x c1f c2f c3f c4f c5f

global @inline function expk{T<:FloatTypes}(d::Double{T})
    q = xrint((d.hi + d.lo)*T(MLN2E))
    s = ddadd2(d, q * -LN2U(T))
    s = ddadd2(s, q * -LN2L(T))
    s = normalize(s)
    u =_expk(s.hi)
    t = ddadd(s, ddmul(ddsqu(s), u))
    t = ddadd(T(1), t)
    ldexpk(t.hi + t.lo, q)
end

global @inline function expk2{T<:FloatTypes}(d::Double{T})
    q = xrint((d.hi + d.lo)*T(MLN2E))
    s = ddadd2(d, q * -LN2U(T))
    s = ddadd2(s, q * -LN2L(T))
    u =_expk(s.hi)
    t = ddadd(s, ddmul(ddsqu(s), u))
    t = ddadd(T(1), t)
    scale(t, pow2i(T, q))
end
end

let
const c8d = 0.134601987501262130076155
const c7d = 0.132248509032032670243288
const c6d = 0.153883458318096079652524
const c5d = 0.181817427573705403298686
const c4d = 0.222222231326187414840781
const c3d = 0.285714285651261412873718
const c2d = 0.400000000000222439910458
const c1d = 0.666666666666666371239645

const c4f = 0.237159967422485351562500f0
const c3f = 0.285279005765914916992188f0
const c2f = 0.400005519390106201171875f0
const c1f = 0.666666567325592041015625f0

global @inline _logk(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d
global @inline _logk(x::Float32) = @horner x c1f c2f c3f c4f

global @inline function logk{T<:FloatTypes}(d::T)
    e  = ilogbp1(d * T(M1SQRT2))
    m  = ldexpk(d,-e)
    x  = dddiv(ddadd2(-T(1), m), ddadd2(T(1), m))
    x2 = ddsqu(x)
    t  =_logk(x2.hi)
    ddadd2(ddmul(MDLN2(T), T(e)), ddadd2(scale(x, T(2.0)), ddmul(ddmul(x2, x), t)))
end

global @inline function logk2{T<:FloatTypes}(d::Double{T})
    e  = ilogbp1(d.hi * T(M1SQRT2))
    m  = scale(d, pow2i(T,-e))
    x  = dddiv(ddadd2(m, -T(1)), ddadd2(m, T(1)))
    x2 = ddsqu(x)
    t  =_logk(x2.hi)
    ddadd2(ddmul(MDLN2(T), T(e)), ddadd2(scale(x, T(2.0)), ddmul(ddmul(x2, x), t)))
end
end
