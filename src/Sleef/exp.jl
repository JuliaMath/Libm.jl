xldexp(x::Float64, q::Int32) = ldexpk(x, q)

function xexp2(a::Float64)
    u = expk(ddmul_d2_d2_d(Double2(0.69314718055994528623, 2.3190468138462995584e-17), a))
    a > 1023   && (u = Inf)
    a === -Inf && (u = 0.0)
  return u
end

function xexp10(a::Float64)
    u = expk(ddmul_d2_d2_d(Double2(2.3025850929940459011, -2.1707562233822493508e-16), a))
    a > 308    && (u = Inf)
    a === -Inf && (u = 0.0)
    return u
end

function xexpm1(a::Float64)
    d = ddadd2_d2_d2_d(expk2(Double2(a, 0.0)), -1.0)
    x = d.x + d.y
    a > 700 && (x = Inf)
    a < -0.36043653389117156089696070315825181539851971360337e+2 && (x = -1.0)
    return x
end

let
global xexp
const c11 = 2.08860621107283687536341e-09
const c10 = 2.51112930892876518610661e-08
const c9  = 2.75573911234900471893338e-07
const c8  = 2.75572362911928827629423e-06
const c7  = 2.4801587159235472998791e-05
const c6  = 0.000198412698960509205564975
const c5  = 0.00138888888889774492207962
const c4  = 0.00833333333331652721664984
const c3  = 0.0416666666666665047591422
const c2  = 0.166666666666666851703837
const c1  = 0.5

function xexp(d::Float64)
    q = xrint(d*R_LN2)
    s = mla(q,-L2U, d)
    s = mla(q,-L2L, s)
    u = @horner s c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11
    u = s * s * u + s + 1
    u = ldexpk(u, q)
    d === -Inf && (u = 0.0)
    return u
end
end
