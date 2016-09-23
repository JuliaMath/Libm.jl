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

## constants

const LN2   = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01 #log(2)
const LOG2E = 1.442695040888963407359924681001892137426645954152985934135449406931 # log2(e)

const M1PI = 0.318309886183790671538 # 1/pi
const M2PI = 0.636619772367581343076 # 2/pi
const MPI  = 3.14159265358979323846  # pi
const M4PI = 1.273239544735162542821171882678754627704620361328125 # 4/pi with round down

## constants (dispatch)

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

# utility functions used by the private math functions in priv.jl

@pure exponent_max{T<:FloatTypes}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))

# _sign emits better native code than sign but does not properly handle the Inf/NaN cases
@inline _sign{T<:FloatTypes}(d::T) =  flipsign(one(T), d) 

@inline xrint{T<:FloatTypes}(x::T) = unsafe_trunc(Int, ifelse(x < 0, x - T(0.5), x + T(0.5)))
# @inline xrintf{T<:FloatTypes}(x::T) = trunc(ifelse(x < 0, x - T(0.5), x + T(0.5)))

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline pow2i{T<:FloatTypes}(::Type{T}, q::Int) = integer2float(T, q + exponent_bias(T))

# sqrt without the domain checks that we don't need since we handle the checks ourselves
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T,x)))

end