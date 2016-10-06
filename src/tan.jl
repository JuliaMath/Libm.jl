# Split 4/pi into four parts (each is 26 bits)
PI4A(::Type{Float64}) = 0.78539816290140151978 
PI4B(::Type{Float64}) = 4.9604678871439933374e-10
PI4C(::Type{Float64}) = 1.1258708853173288931e-18
PI4D(::Type{Float64}) = 1.7607799325916000908e-27

PI4A(::Type{Float32}) = 0.78515625f0
PI4B(::Type{Float32}) = 0.00024187564849853515625f0
PI4C(::Type{Float32}) = 3.7747668102383613586f-08
PI4D(::Type{Float32}) = 1.2816720341285448015f-12

let
const c15d =  1.01419718511083373224408e-05
const c14d = -2.59519791585924697698614e-05
const c13d =  5.23388081915899855325186e-05
const c12d = -3.05033014433946488225616e-05
const c11d =  7.14707504084242744267497e-05
const c10d =  8.09674518280159187045078e-05
const c9d  =  0.000244884931879331847054404
const c8d  =  0.000588505168743587154904506
const c7d  =  0.001456127889228124279788480
const c6d  =  0.003592087438369066191429240
const c5d  =  0.008863239443624016181133560
const c4d  =  0.021869488285384638959207800
const c3d  =  0.053968253978129841763600200
const c2d  =  0.133333333333125941821962000
const c1d  =  0.333333333333334980164153000

const c7f =  0.00446636462584137916564941
const c6f = -8.3920182078145444393158e-05
const c5f =  0.0109639242291450500488281
const c4f =  0.0212360303848981857299805
const c3f =  0.0540687143802642822265625
const c2f =  0.133325666189193725585938
const c1f =  0.33333361148834228515625

global @inline _tan(x) = @horner_split_oftype x c1d c2d c3d c4d c5d c6d c7d c8d c9d c10d c11d c12d c13d c14d c15d
global @inline _tan(x::SFloat) = @horner_split_oftype x c1f c2f c3f c4f c5f c6f c7f

# does not handle -0.0
global function tan{T}(dd::T)
    d = abs(dd)
    q = round(d*T(M2PI))
    x = muladd(q, -PI4A(T)*2, d)
    x = muladd(q, -PI4B(T)*2, x)
    x = muladd(q, -PI4C(T)*2, x)
    x = muladd(q, -PI4D(T)*2, x)
    n = _trunc(q)
    x *= T(1.0) - T(2.0)*(n & 1) # n & 1 != 0 && (x = -x)
    
    s = x*x
    u =_tan(s)
    u = muladd(s,u*x,x)

     #n & 1 != 0 && (u = 1/u), the next line does the same with no branch
     u = (1-(n & 1))*u + (n & 1)*(T(1.0)/(u+realmin(T))) # branchless version, add eps to prevent 1/0 exceptional case
    return flipsign(u,dd)
end
end
