module Sleef

export xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
    xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp

# Higher accuracy functions
export xatan2_u1, xasin_u1, xacos_u1, xatan_u1, xsin_u1, xcos_u1, xsincos_u1, xtan_u1, xcbrt_u1, xlog_u1

# Alias for supported floating point types
typealias FloatTypes Union{Float32,Float64}

# Split 4/pi into four parts (each is 26 bits)
const PI4A = 0.78539816290140151978 
const PI4B = 4.9604678871439933374e-10
const PI4C = 1.1258708853173288931e-18
const PI4D = 1.7607799325916000908e-27
const M4PI = 1.273239544735162542821171882678754627704620361328125 # 4/pi with round down

# Split log(2) into upper and lower parts
const LN2U = 0.69314718055966295651160180568695068359375
const LN2L = 0.28235290563031577122588448175013436025525412068e-12

const LN2   = 6.931471805599453094172321214581765680755001343602552541206800094933936219696955e-01
const LOG2E = 1.442695040888963407359924681001892137426645954152985934135449406931 # log2(e)

const M1PI = 0.318309886183790671538 # 1/pi
const M2PI = 0.636619772367581343076 # 2/pi
const MPI  = 3.14159265358979323846  # pi

using Base: Math.@horner, significand_bits, exponent_bits, exponent_bias, exponent_mask, @pure

include("Sleef/double.jl")
include("Sleef/priv.jl") # private math functions

# exported math functions
include("Sleef/exp.jl")  # exponential functions
include("Sleef/log.jl")  # logarithmic functions
include("Sleef/trig.jl") # trigonometric and inverse trigonometric functions
include("Sleef/hyp.jl")  # hyperbolic and inverse hyperbolic functions
include("Sleef/misc.jl") # miscallenous math functions including pow and cbrt

# utility functions used by the private math functions in priv.jl

@pure exponent_max{T<:AbstractFloat}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))

@inline sign{T<:FloatTypes}(d::T) =  copysign(one(T), d) # emits better native code than Base.sign

@inline xrint{T<:FloatTypes}(x::T) = x < 0 ? unsafe_trunc(Int, x - T(0.5)) : unsafe_trunc(Int, x + T(0.5))

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline pow2i(::Type{Float64}, q::Int) = integer2float(Float64, q + exponent_bias(Float64))
@inline pow2i(::Type{Float32}, q::Int) = integer2float(Float32, q + exponent_bias(Float64))

@inline upper(d::Float64) = reinterpret(Float64, reinterpret(UInt64, d) & 0xfffffffff8000000) # clear lower 27 bits (leave upper 26 bits)
@inline upper(d::Float32) = reinterpret(Float32, reinterpret(UInt32, d) & 0xfffff000) # clear lowest 12 bits (leave upper 12 bits)

# sqrt without the domain checks that we don't need since we handle the checks ourselves
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T,x)))

end