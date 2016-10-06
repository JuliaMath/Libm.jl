module Amal

export exp, exp2, atan, tan

using Base.Math: @horner

typealias Float Union{Float16,Float32,Float64} # concrete float types
typealias HFloat Union{Float64} # large floats
typealias SFloat Union{Float16,Float32} # small floats

# Float128
const PI = 3.141592653589793238462643383279502884197169
const PI2 = 1.570796326794896619231321691639751442098585
const PI4 =  0.7853981633974483096156608458198757210492923
const M2PI = 0.63661977236758138243288840385503135621547698974609375
const LOGE2 =  0.6931471805599453094172321214581765680755001
const LOG2E =  1.442695040888963407359924681001892137426646


MAXLOG(::Type{Float64}) =  7.09782712893383996732e2   # log 2^1023*(2-2^-52)
MINLOG(::Type{Float64}) = -7.451332191019412076235e2  # log 2^-1075 (one less than the min exponent, since we can actually sqeeze just a bit more from the exp function)

MAXLOG(::Type{Float32}) =  88.72283905206835f0      #  log 2^127*(2-2^-23)
MINLOG(::Type{Float32}) = -103.972077083991796412584818218726485211325f0 # log 2^-149 (one less than min exp)

MAXLOG(::Type{Float16}) = Float16(11.09)  # log 2^15*(2-2^-10)
MINLOG(::Type{Float16}) = Float16(-17.32868) # log 2^-25 (one less than min exp)



# Similar to @horner, but converts coefficients to same type as x
macro horner_oftype(x, p...)
    val = gensym()
    ex = :(oftype($x,$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd($val, $ex, oftype($x,$(esc(p[i])))))
    end
    Expr(:block, :($val = $(esc(x))), ex)
end

# Similar to @horner, but split into even and odd coefficients.
macro horner_split_oftype(x,p...)
    t1 = gensym()
    t2 = gensym()
    blk = quote
        $t1 = $(esc(x))
        $t2 = $(esc(x)) * $(esc(x))
    end
    n = length(p)
    p0 = :(oftype($x,$(esc(p[1]))))
    if isodd(n)
        ex_o = :(oftype($x,$(esc(p[end-1]))))
        ex_e = :(oftype($x,$(esc(p[end]))))
        for i = n-3:-2:2
            ex_o = :(muladd($(t2), $ex_o, oftype($x,$(esc(p[i])))))
        end
        for i = n-2:-2:2
            ex_e = :(muladd($(t2), $ex_e, oftype($x,$(esc(p[i])))))
        end
    elseif iseven(n)
        ex_o = :(oftype($x,$(esc(p[end]))))
        ex_e = :(oftype($x,$(esc(p[end-1]))))
        for i = n-2:-2:2
            ex_o = :(muladd($(t2), $ex_o, oftype($x,$(esc(p[i])))))
        end
        for i = n-3:-2:2
            ex_e = :(muladd($(t2), $ex_e, oftype($x,$(esc(p[i])))))
        end
    end
    push!(blk.args,:($(p0) + $(t1)*$(ex_o) + $(t2)*$(ex_e)))
    blk
end


# Float16 hacks
import Base: unsafe_trunc
unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x))

#∧{T}(x::Real, ::Type{T}) = T(x) #or \leftthreetimes


using Base: significand_bits, exponent_bias, exponent_mask, exponent_one
import Base: exponent_bias

_trunc{T<:Float}(x::T) = unsafe_trunc(Ti(T),x)

Ti(::Type{Float64}) = Int64
Ti(::Type{Float32}) = Int32
Ti(::Type{Float16}) = Int16
Ti(::Type{Float64},x) = Int64(x)
Ti(::Type{Float32},x) = Int32(x)
Ti(::Type{Float16},x) = Int16(x)

integer2float(::Type{Float64}, m::Int64) = reinterpret(Float64, (m % Int64) << significand_bits(Float64))
integer2float(::Type{Float32}, m::Int32) = reinterpret(Float32, (m % Int32) << significand_bits(Float32))
integer2float(::Type{Float16}, m::Int16) = reinterpret(Float16, (m % Int16) << significand_bits(Float16))
float2integer(d::Float64) = (reinterpret(Int64, d) >> significand_bits(Float64))
float2integer(d::Float32) = (reinterpret(Int32, d) >> significand_bits(Float32))
float2integer(d::Float16) = (reinterpret(Int16, d) >> significand_bits(Float16))

exponent_bias{T<:Float}(::Type{T}) = Ti(T, exponent_one(T) >> significand_bits(T))
exponent_max{T<:Float}(::Type{T}) = Ti(T, exponent_mask(T) >> significand_bits(T))


@inline split_exponent(::Type{Float64}, q::Int64) = _split_exponent(q, UInt64(9), UInt64(31), UInt64(2))
@inline split_exponent(::Type{Float32}, q::Int32) = _split_exponent(q, UInt32(6), UInt32(31), UInt32(2))
@inline split_exponent(::Type{Float16}, q::Int16) = _split_exponent(q, UInt16(3), UInt16(31), UInt16(2))
@inline function _split_exponent(q, n, v, offset)
    m = q >> v
    m = (((m + q) >> n) - m) << (n-offset)
    q = q - (m << offset)
    m, q
end

"""
    ldexpk(a::Float, n::Int) -> Float

Computes `a × 2^n`.
"""
@inline function ldexpk{T<:Float}(x::T, q::Integer)
    bias = exponent_bias(T)
    emax = exponent_max(T)
    m, q = split_exponent(T,q)
    m += bias
    m = ifelse(m < 0, Ti(T,0), m)
    m = ifelse(m > emax, emax, m)
    q += bias
    u = integer2float(T, m)
    x = x*u*u*u*u
    u = integer2float(T, q)
    x*u
end


# Float128
# C1 + C2 = -ln 2 
C1{T<:Float}(::Type{T}) = T(6.93145751953125e-1)
C2{T<:Float}(::Type{T}) = T(1.428606820309417232121458176568075500134e-6)

C1{T<:SFloat}(::Type{T}) =  T(0.693359375)
C2{T<:SFloat}(::Type{T}) =  T(-2.12194440e-4)


include("exp.jl")
include("exp2.jl")
include("atan.jl")
include("tan.jl")

end