xldexp{T<:FloatTypes}(x::T, q::Int) = ldexpk(x, q)

function xexp2(a::Float64)
    u = expk(ddmul(Double(0.69314718055994528623, 2.3190468138462995584e-17), a))
    a > 1023   && (u = Inf)
    a === -Inf && (u = 0.0)
  return u
end

function xexp10(a::Float64)
    u = expk(ddmul(Double(2.3025850929940459011, -2.1707562233822493508e-16), a))
    a > 308    && (u = Inf)
    a === -Inf && (u = 0.0)
    return u
end

function xexpm1(a::Float64)
    d = ddadd2(expk2(Double(a)), -1.0)
    x = d.hi + d.lo
    a > 700 && (x = Inf)
    a < -0.36043653389117156089696070315825181539851971360337e+2 && (x = -1.0)
    return x
end

let
global xexp
const c11d = 2.08860621107283687536341e-09
const c10d = 2.51112930892876518610661e-08
const c9d  = 2.75573911234900471893338e-07
const c8d  = 2.75572362911928827629423e-06
const c7d  = 2.4801587159235472998791e-05
const c6d  = 0.000198412698960509205564975
const c5d  = 0.00138888888889774492207962
const c4d  = 0.00833333333331652721664984
const c3d  = 0.0416666666666665047591422
const c2d  = 0.166666666666666851703837
const c1d  = 0.5

const c5f = 0.00136324646882712841033936f0
const c4f = 0.00836596917361021041870117f0
const c3f = 0.0416710823774337768554688f0
const c2f = 0.166665524244308471679688f0
const c1f = 0.499999850988388061523438f0

@inline _xexp(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d  
# @inline _xexp(x::Float32) = @horner x c1f c2f c3f c4f c5f

function xexp{T<:Float64}(d::T)
    q = xrint(d*T(LOG2E))
    s = muladd(q,-LN2U, d)
    s = muladd(q,-LN2L, s)
    u = _xexp(s)
    u = s*s*u + s + 1
    u = ldexpk(u, q)
    d === -Inf && (u = 0.0)
    return u
end
end
