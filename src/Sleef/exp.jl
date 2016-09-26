xldexp(x::FloatTypes, q::Int) = ldexpk(x, q)

function xexp2{T<:FloatTypes}(a::T)
    u = expk(ddmul(MDLN2(T), a))
    ispinf(a) && (u = typemax(T))
    isninf(a) && (u = T(0))
    return u
end

function xexp10{T<:FloatTypes}(a::T)
    u = expk(ddmul(MDLN10(T), a))
    ispinf(a) && (u = typemax(T))
    isninf(a) && (u = T(0))
    return u
end

# over/under flow argument values
_over_xexpm1(::Type{Float64}) =  700.0
_undr_xexpm1(::Type{Float64}) = -0.36043653389117156089696070315825181539851971360337e2
_over_xexpm1(::Type{Float32}) =  88f0
_undr_xexpm1(::Type{Float32}) = -0.15942385152878742116596338793538061065739925620174f2

function xexpm1{T<:FloatTypes}(a::T)
    d = ddadd2(expk2(Double(a)), -T(1))
    x = d.hi + d.lo
    a > _over_xexpm1(T) && (x =  typemax(T))
    a < _undr_xexpm1(T) && (x = -T(1))
    return x
end

let
global xexp
const c11d = 2.08860621107283687536341e-09
const c10d = 2.51112930892876518610661e-08
const c9d  = 2.75573911234900471893338e-07
const c8d  = 2.75572362911928827629423e-06
const c7d  = 2.48015871592354729987910e-05
const c6d  = 0.000198412698960509205564975
const c5d  = 0.001388888888897744922079620
const c4d  = 0.008333333333316527216649840
const c3d  = 0.041666666666666504759142200
const c2d  = 0.166666666666666851703837000
const c1d  = 0.500000000000000000000000000

const c5f = 0.00136324646882712841033936f0
const c4f = 0.00836596917361021041870117f0
const c3f = 0.04167108237743377685546880f0
const c2f = 0.16666552424430847167968800f0
const c1f = 0.49999985098838806152343800f0

global @inline _xexp(x::Float64) = @horner x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d  
global @inline _xexp(x::Float32) = @horner x c1f c2f c3f c4f c5f

function xexp{T<:FloatTypes}(d::T)
    q = xrint(d*T(MLN2E))
    s = muladd(q, -LN2U(T), d)
    s = muladd(q, -LN2L(T), s)
    u =_xexp(s)
    u = s*s*u + s + T(1)
    u = ldexpk(u, q)
    isninf(d) && (u = T(0))
    return u
end
end
