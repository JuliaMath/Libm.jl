let
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
const c1d  = 0.50

const c5f = 0.00136324646882712841033936
const c4f = 0.00836596917361021041870117
const c3f = 0.0416710823774337768554688
const c2f = 0.166665524244308471679688
const c1f = 0.499999850988388061523438

global @inline _exp2(x) = @horner_oftype x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d  
global @inline _exp2(x::SFloat) = @horner_oftype x c1f c2f c3f c4f c5f 

global function exp2{T}(x::T)
    q = _trunc(round(x)) # truncation will give us automatic Inf handling
    t::T = x - q  # for float16 bug on 0.5
    s = t*C1(T)
    s = muladd(t,C2(T), s)

    u =_exp2(s)
    u = muladd(s,s*u,s) + T(1.0)
    u = ldexpk(u,q)
    
    x == T(-Inf) && (u = T(0.0))
    return u
end
end

