module Cephes

using Base.Math: @horner

typealias Float Union{Float16,Float32,Float64} # concrete float types only please
typealias HFloat Union{Float64} # small floats
typealias SFloat Union{Float16,Float32} # small floats

# Float128
const PI = 3.141592653589793238462643383279502884197169
const PI2 = 1.570796326794896619231321691639751442098585
const PI4 =  0.7853981633974483096156608458198757210492923
const LOGE2 =  0.6931471805599453094172321214581765680755001
const LOG2E =  1.442695040888963407359924681001892137426646


MAXLOG(::Type{Float64}) =  7.09782712893383996732e2   # log 2^1023*(2-2^-52)
MINLOG(::Type{Float64}) = -7.451332191019412076235e2  # log 2^-1075 (one less than min exp)

MAXLOG(::Type{Float32}) =  88.72283905206835f0      #  log 2^127*(2-2^-23)
MINLOG(::Type{Float32}) = -103.972077083991796412584818218726485211325f0 # log 2^-149 (one less than min exp)

MAXLOG(::Type{Float16}) = Float16(11.09)  # log 2^15*(2-2^-10)
MINLOG(::Type{Float16}) = Float16(-17.32868) # log 2^-25 (one less than min exp)

# float16 hacks
import Base: unsafe_trunc
unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x))

# Similar to @horner, but converts coefficients to same type as x
macro horner_oftype(x, p...)
    val = gensym()
    ex = :(oftype($x,$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd($val, $ex, oftype($x,$(esc(p[i])))))
    end
    Expr(:block, :($val = $(esc(x))), ex)
end

#∧{T}(x::Real, ::Type{T}) = T(x) #or \leftthreetimes

include("atan.jl")
include("exp.jl")

end