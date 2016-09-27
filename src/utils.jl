## utility functions mainly used by the private math functions in priv.jl

function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)),nextfloat(one(T)),-nextfloat(one(T),2)) != zero(T))
end
is_fma_fast() = is_fma_fast(Float64) && is_fma_fast(Float32)


@inline exponent_max{T<:FloatTypes}(::Type{T}) = Int(exponent_mask(T) >> significand_bits(T))

# _sign emits better native code than sign but does not properly handle the Inf/NaN cases
@inline _sign{T<:FloatTypes}(d::T) =  flipsign(T(1), d)

@inline xrint{T<:FloatTypes}(x::T) = unsafe_trunc(Int, ifelse(x < 0, x - T(0.5), x + T(0.5)))
# @inline xrintf{T<:FloatTypes}(x::T) = trunc(ifelse(x < 0, x - T(0.5), x + T(0.5)))

@inline integer2float(::Type{Float64}, m::Int) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
@inline integer2float(::Type{Float32}, m::Int) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))

@inline float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64)) % Int
@inline float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32)) % Int

@inline pow2i{T<:FloatTypes}(::Type{T}, q::Int) = integer2float(T, q + exponent_bias(T))

# sqrt without the domain checks which we don't need since we handle the checks ourselves
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T,x)))

@inline ispinf{T<:FloatTypes}(x::T) = x == typemax(T)
@inline isninf{T<:FloatTypes}(x::T) = x == typemin(T)
