# This file includes the exported function xilogb, xlog_u1, xlog1p, and xlog

"""
    xilogb(x::FloatTypes) -> FloatTypes

Returns the integral part of the logarithm of `abs(x)`, using base 2 for the logarithm; in other
words this returns the binary exponent of `x` so that

    x = significand \times 2^exponent,

where `significand \in [1, 2)`

* Exceptional cases (where `Int` is the wordsize, which is either `Int64` or `Int32`)
    * `x = 0`    returns `typemin(Int)`
    * `x = Inf`  returns `typemax(Int)`
    * `x = NaN`  returns `typemax(Int)`
"""
function xilogb{T<:FloatTypes}(d::T)
    e = ilogbp1(abs(d)) - 1
    e = d == 0   ? typemin(Int) : e
    e = isinf(d) ? typemax(Int) : e
    e = isnan(d) ? typemax(Int) : e
    return e 
end

function xlog_u1{T<:FloatTypes}(d::T)
    s = logk(d)
    x = s.hi + s.lo
    isinf(d) && (x =  T(Inf))
    d < 0    && (x =  T(NaN))
    d == 0   && (x = -T(Inf))
    return x
end


_xlog10_c0(::Type{Float32}) = Double(0.43429449200630187988f0, -1.0103050118726031315f-08)
_xlog10_c0(::Type{Float64}) = Double(0.43429448190325176116, 6.6494347733425473126e-17)

function xlog10{T<:FloatTypes}(a::T)
    d = ddmul(logk(a), _xlog10_c0(T))
    x = d.hi + d.lo
    isinf(a) && (x =  T(Inf))
    a < 0    && (x =  T(NaN))
    a == 0   && (x = -T(Inf))
    return x
end


function xlog1p{T<:FloatTypes}(a::T)
    d = logk2(ddadd2(a, T(1)))
    x = d.hi + d.lo
    isinf(a) && (x =  T(Inf))
    a < -1   && (x =  T(NaN))
    a == -1  && (x = -T(Inf))
    return x
end




# First we split the argument to its mantissa `m` and integer exponent `e` so that `d = m \times 2^e`,
# where `m \in [0.5, 1)` then we apply the polynomial approximant on this reduced argument `m` before
# putting back the exponent in. This first part is done with the help of the private function
# `ilogbp1(x)` and we put the exponent back using

#     `\log(m \times 2^e) = \log(m) + \log 2^e =  \log(m) + e\times LN2

# The polynomial we evaluate is based on coefficients from

#     `log_2(x) = 2\sum_{n=0}^\infty \frac{1}{2n+1} \bigl(\frac{x-1}{x+1}^{2n+1}\bigr)`

# That being said, since this converges faster when the argument is close to 1, we multiply  `m` by
# `2` and subtract 1 for the exponent `e` when `m` is less than `sqrt(2)/2`

let
global xlog
const c8d = 0.148197055177935105296783
const c7d = 0.153108178020442575739679
const c6d = 0.181837339521549679055568
const c5d = 0.22222194152736701733275
const c4d = 0.285714288030134544449368
const c3d = 0.399999999989941956712869
const c2d = 0.666666666666685503450651
const c1d = 2

const c5f = 0.2371599674224853515625f0
const c4f = 0.285279005765914916992188f0
const c3f = 0.400005519390106201171875f0
const c2f = 0.666666567325592041015625f0
const c1f = 2.0f0

global @inline _xlog(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d
global @inline _xlog(x::Float32) = @horner x c1f c2f c3f c4f c5f

function xlog{T<:FloatTypes}(d::T)
    e = ilogbp1(d*T(M1SQRT2))
    m = ldexpk(d,-e)
    x = (m-1)/(m+1)
    x2 = x*x
    t = _xlog(x2)
    x = muladd(x, t, T(LN2)*e)
    isinf(d) && (x =  T(Inf))
    d < 0    && (x =  T(NaN))
    d == 0   && (x = -T(Inf))
    return x
end
end
