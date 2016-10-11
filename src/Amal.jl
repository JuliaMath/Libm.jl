module Amal

using Base: significand_bits, exponent_bias, exponent_mask, exponent_one

typealias FloatTypes Union{Float16,Float32,Float64}
typealias IntTypes Union{Int16,Int32,Int64}

typealias LargeFloatTypes Union{Float64}
typealias SmallFloatTypes Union{Float16,Float32}

inttype(::Type{Float64}) = Int64
inttype(::Type{Float32}) = Int32
inttype(::Type{Float16}) = Int16

floattype(::Type{Int64}) = Float64
floattype(::Type{Int32}) = Float32
floattype(::Type{Int16}) = Float16

# convert float to corresponding integer type of the same size
asint{T<:FloatTypes}(::Type{T}, x) = inttype(T)(x)

# reinterpret an integer to the corresponding float of the same size
intasfloat{T<:IntTypes}(m::T) = reinterpret(floattype(T), m << significand_bits(floattype(T)))
# reinterpret a float to the corresponding int of the same size
floatasint{T<:FloatTypes}(d::T) = reinterpret(inttype(T), d) >> significand_bits(T)

exponent_max{T<:FloatTypes}(::Type{T}) = asint(T, exponent_mask(T) >> significand_bits(T))

# import Base: unsafe_trunc
# unsafe_trunc{T<:Integer}(::Type{T}, x::Float16) = trunc(T, Float32(x)) # Float16 hack for v0.5
# unsafe truncate x and return integer of the same size as x
_trunc{T<:FloatTypes}(x::T) = unsafe_trunc(inttype(T), x)



# evaluate p[1] + x * (p[2] + x * (....)), i.e. a polynomial via Horner's rule
# and convert coefficients to same type as x
macro horner_oftype(x, p...)
    @gensym val
    ex = :(oftype($x,$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd($val, $ex, oftype($x,$(esc(p[i])))))
    end
   return Expr(:block, :($val = $(esc(x))), ex)
end

function _numeric(T,ex) # T is a Symbol
    isa(ex, Symbol) && return ex
    if isa(ex, Expr)
        ex.args = map(x -> isa(x, Number) ? :($T($x)) : _numeric(T,x), ex.args)
        return ex
    end
end

macro oftype(ex)
    if is(ex.head,:(=)) && is(ex.args[1].head,:call) || is(ex.head,:function) 
        if isa(ex.args[1].args[1],Symbol) && is(ex.args[1].args[2].head,:(::))
            type_param = gensym("T")
            typ = ex.args[1].args[2].args[2]
            fun_name = ex.args[1].args[1]
            ex.args[1].args[1] = Expr(:curly, fun_name, Expr(:<:, type_param ,typ))
            ex.args[1].args[2].args[2] = type_param
            ex.args[2]  = _numeric(type_param, ex.args[2])
            return esc(ex)
        elseif is(ex.args[1].args[1].head,:curly) && isa(ex.args[1].args[1].args[2],Expr)
            type_param = ex.args[1].args[1].args[2].args[1]
            ex.args[2]  = _numeric(type_param, ex.args[2])
            return esc(ex)
        elseif is(ex.args[1].args[1].head,:curly) 
            type_param = ex.args[1].args[1].args[2]
            ex.args[2]  = _numeric(type_param, ex.args[2])
            return esc(ex)
        end
    end
end


# constants

const LOG2E =  1.442695040888963407359924681001892137426646

LN2U{T<:LargeFloatTypes}(::Type{T}) = T(0.69314718055966295651160180568695068359375)
LN2L{T<:LargeFloatTypes}(::Type{T}) = T(0.28235290563031577122588448175013436025525412068e-12)

LN2U{T<:SmallFloatTypes}(::Type{T}) = T(0.693145751953125)
LN2L{T<:SmallFloatTypes}(::Type{T}) = T(1.428606765330187045e-06)


include("ldexp.jl")

include("exp.jl") # exp, exp2
# include("poly/exp.jl");
# include("poly/exp2.jl")


end