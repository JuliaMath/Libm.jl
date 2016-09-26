module Sleef

export xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
    xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp
# Higher accuracy functions
export xatan2_u1, xasin_u1, xacos_u1, xatan_u1, xsin_u1, xcos_u1, xsincos_u1, xtan_u1, xcbrt_u1, xlog_u1

# Alias for supported floating point types
typealias FloatTypes Union{Float32,Float64}

using Base: Math.@horner, significand_bits, exponent_bits, exponent_bias, exponent_mask, @pure

function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)),nextfloat(one(T)),-nextfloat(one(T),2)) != zero(T))
end
is_fma_fast() = is_fma_fast(Float64) && is_fma_fast(Float32)

## constants (refactor lator)

const MLN2  = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01 # log(2)
const MLN2E = 1.442695040888963407359924681001892137426645954152985934135449406931109219181187     # log2(e)

const MPI  = 3.141592653589793238462643383279502884197169399375105820974944592307816406286198     # pi
const MPI2 = 1.570796326794896619231321691639751442098584699687552910487472296153908203143099     # pi/2
const MPI4 = 7.853981633974483096156608458198757210492923498437764552437361480769541015715495e-01 # pi/4
const M1PI = 3.183098861837906715377675267450287240689192914809128974953346881177935952684543e-01 # 1/pi
const M2PI = 6.366197723675813430755350534900574481378385829618257949906693762355871905369086e-01 # 2/pi
const M4PI = 1.273239544735162686151070106980114896275677165923651589981338752471174381073817     # 4/pi

const M1SQRT2 = 7.07106781186547524400844362104849039284835937688474036588339868995366239231051e-01 # 1/sqrt(2)

## constants (dispatch)

MDLN10E(::Type{Float64}) = Double(0.4342944819032518, 1.098319650216765e-17) # log10(e)
MDLN10E(::Type{Float32}) = Double(0.4342945f0, -1.010305f-8)

MDLN10(::Type{Float64}) = Double(2.302585092994046, -2.1707562233822494e-16) # log(10)
MDLN10(::Type{Float32}) = Double(2.3025851f0, -3.1975436f-8)

MDLN2(::Type{Float64}) = Double(0.6931471805599453, 2.3190468138462996e-17)  # log(2)
MDLN2(::Type{Float32}) = Double(0.6931472f0, -1.9046542f-9)

MDPI(::Type{Float64})  = Double(3.141592653589793, 1.2246467991473532e-16) # pi
MDPI(::Type{Float32})  = Double(3.1415927f0, -8.742278f-8)
MDPI2(::Type{Float64}) = Double(1.5707963267948966, 6.123233995736766e-17) # pi/2
MDPI2(::Type{Float32}) = Double(1.5707964f0, -4.371139f-8)

# Split 4/pi into four parts (each is 26 bits)
PI4A(::Type{Float64}) = 0.78539816290140151978 
PI4B(::Type{Float64}) = 4.9604678871439933374e-10
PI4C(::Type{Float64}) = 1.1258708853173288931e-18
PI4D(::Type{Float64}) = 1.7607799325916000908e-27

PI4A(::Type{Float32}) = 0.78515625f0
PI4B(::Type{Float32}) = 0.00024187564849853515625f0
PI4C(::Type{Float32}) = 3.7747668102383613586f-08
PI4D(::Type{Float32}) = 1.2816720341285448015f-12

# Split log(2) into upper and lower parts
LN2U(::Type{Float64}) = 0.69314718055966295651160180568695068359375
LN2L(::Type{Float64}) = 0.28235290563031577122588448175013436025525412068e-12

LN2U(::Type{Float32}) = 0.693145751953125f0
LN2L(::Type{Float32}) = 1.428606765330187045f-06


include("Sleef/double.jl") # Dekker style Double implementation
include("Sleef/priv.jl")   # private math functions
include("Sleef/exp.jl")    # exponential functions
include("Sleef/log.jl")    # logarithmic functions
include("Sleef/trig.jl")   # trigonometric and inverse trigonometric functions
include("Sleef/hyp.jl")    # hyperbolic and inverse hyperbolic functions
include("Sleef/misc.jl")   # miscallenous math functions including pow and cbrt

## utility functions mainly used by the private math functions in priv.jl

@inline exponent_max{T<:FloatTypes}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))

# _sign emits better native code than sign but does not properly handle the Inf/NaN cases
@inline _sign{T<:FloatTypes}(d::T) =  flipsign(one(T), d) 

@inline xrint{T<:FloatTypes}(x::T) = unsafe_trunc(Int, ifelse(x < 0, x - T(0.5), x + T(0.5)))
# @inline xrintf{T<:FloatTypes}(x::T) = trunc(ifelse(x < 0, x - T(0.5), x + T(0.5)))

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline pow2i{T<:FloatTypes}(::Type{T}, q::Int) = integer2float(T, q + exponent_bias(T))

# sqrt without the domain checks which we don't need since we handle the checks ourselves
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T,x)))

@inline ispinf{T<:FloatTypes}(x::T) = x ==  T(Inf)
@inline isninf{T<:FloatTypes}(x::T) = x == -T(Inf)


end