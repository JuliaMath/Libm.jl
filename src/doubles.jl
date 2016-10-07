# module Doubles
# export Double, twosum_fast, twosum, twoprod, scale, abs
export Double
import Base: +, -, *, /, sqrt, abs, scale, ==, convert, promote_rule, show

Float = Union{Float16,Float32,Float64}
abstract AbstractDouble{T} <: Real
immutable Single{T<:Float} <: AbstractDouble{T}
    hi::T
end
immutable Double{T<:Float} <: AbstractDouble{T}
    hi::T
    lo::T
end
Double(x) = Single(x)
Double{T}(x::Tuple{T,T}) = Double(x...)
DoubleNorm(x,y) = Double(twosum_fast(x,y)) # normalize the double

hibits{T}(x::Double{T}) = x.hi

### promotions and conversions ###

convert{T<:Float}(::Type{T}, x::Single) = x.hi
convert{T<:Float}(::Type{T}, x::Double) = x.hi + x.lo
convert(::Type{BigFloat}, x::Single) = big(x.hi)
convert(::Type{BigFloat}, x::Double) = big(x.hi) + big(x.lo)

# The following hack promotes the float types to AbstractDouble so that 
# float types get properly converted to a Single type.
# We need this since we do not want to promote floats to Double since
# we want to dispatch to methods with the Single type for effiency reasons.
# Similar for the other types.
convert{T}(::Type{AbstractDouble{T}}, z::T)               = Single(z)
convert{T}(::Type{AbstractDouble{T}}, z::Type{Single{T}}) = z
convert{T}(::Type{AbstractDouble{T}}, z::Type{Double{T}}) = z

convert{T}(::Type{AbstractDouble{T}}, z::Real) = Single(convert(T,z)) # fall back conversion


promote_rule{T}(::Type{Single{T}}, ::Type{T})         = AbstractDouble{T}
promote_rule{T}(::Type{Double{T}}, ::Type{T})         = AbstractDouble{T}
promote_rule{T}(::Type{Double{T}}, ::Type{Single{T}}) = AbstractDouble{T}

promote_rule{T,R<:Real}(::Type{Single{T}}, ::Type{R}) = AbstractDouble{T}  # fall back conversion
promote_rule{T,R<:Real}(::Type{Double{T}}, ::Type{R}) = AbstractDouble{T}  # fall back conversion


### utility functions

# the following are only used for non fma systems
@inline trunclo(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xffff_ffff_f800_0000) # clear lower 27 bits (leave upper 26 bits)
@inline trunclo(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xffff_f000) # clear lowest 12 bits (leave upper 12 bits)
@inline function splitprec(x::Float)
    hx = trunclo(x)
    hx, x-hx
end


### basic error free float arithmetic ###

# fast two-sum addition, if |x| ≥ |y|
@inline function twosum_fast{T<:Float}(x::T, y::T)
    r = x + y
    r, y + (x - r)
end

# two-sum addition
@inline function twosum{T<:Float}(x::T, y::T)
    r = x + y
    v = r - x
    r, (y - v) + (x - (r - v))
end

# two-product fma
@inline function _twoprod_fma{T<:Float}(x::T, y::T)
    r = x*y
    r, fma(x,y,-r)
end

# two-product non-fma
@inline function _twoprod{T<:Float}(x::T, y::T)
    hx, lx = splitprec(x)
    hy, ly = splitprec(y)
    z = x*y
    z, ((hx*hy-z) + lx*hy + hx*ly) + lx*ly
end


### Double arithmetic ###

## negation

@inline -(x::Double)  = Double(-x.hi, -x.lo)
@inline -(x::Single)  = Single(-x.hi)


## addition

@inline function +{T}(x::Double{T}, y::Double{T})
    r, s = twosum(x.hi, y.hi)
    Double(r, s + x.lo + y.lo)
end

@inline function +{T}(x::Double{T}, y::Single{T})
    r, s = twosum(x.hi, y.hi)
    Double(r, s + x.lo)
end
@inline +{T}(x::Single{T}, y::Double{T}{T}) = y + x

@inline +{T}(x::Single{T}, y::Single{T}) = Double(twosum(x.hi, y.hi))

## subtraction

@inline function -{T}(x::Double{T}, y::Double{T})
    r, s = twosum(x.hi, -y.hi)
    Double(r, s + (x.lo - y.lo))
end

@inline function -{T}(x::Double{T}, y::Single{T})
    r, s = twosum(x.hi, -y.hi)
    Double(r, s + x.lo)
end

@inline function -{T}(x::Single{T}, y::Double{T})
    r, s = twosum(x.hi, -y.hi)
    Double(r, s + y.lo)
end

@inline -{T}(x::Single{T}, y::Single{T}) = Double(twosum(x.hi, -y.hi))


## multiplication

@inline function *{T}(x::Double{T}, y::Double{T})
    r, s = twoprod(x.hi, y.hi)
    Double(r, s + x.hi*y.lo + x.lo*y.hi)
end

@inline function *{T}(x::Double{T}, y::Single{T})
    z0, z1 = twoprod(x.hi, y.hi)
    Double(z0, z1 + x.lo)
end
@inline *{T}(x::Single{T}, y::Double{T}) = y*x

@inline *{T}(x::Single{T}, y::Single{T}) = Double(twoprod(x.hi, y.hi))


## division

# private two-product div fma helper
@inline function _pdiv_fma{T<:Float}(x::T, y::T)
    ry = 1/y
    r = x*ry
    r, fma(-r,y,x), ry
end

# private two-product div helper
@inline function _pdiv{T<:Float}(x::T, y::T)
    ry = 1/y
    r = x*ry
    hx, lx = splitprec(r)
    hy, ly = splitprec(y)
    r, ((-hx*hy+r*y) - lx*hy - hx*ly) - lx*ly, ry
end

@inline function /{T}(x::Double{T}, y::Double{T})
    r, s, ry = pdiv(x.hi, y.hi)
    Double(r, (s + muladd(-r, y.lo, x.lo))*ry)
end

@inline function /{T}(x::Double{T}, y::Single{T})
    r, s, ry = pdiv(x.hi, y.hi)
    Double(r, (s - x.lo)*ry)
end

@inline function /{T}(x::Single{T}, y::Double{T})
    r, s, ry = pdiv(x.hi, y.hi)
    Double(r, muladd(-r, y.lo, s)*ry)
end

@inline function /{T}(x::Single{T}, y::Single{T})
    r, s, ry = pdiv(x.hi, y.hi)
    Double(r, s*ry)
end

## square root

# fast sqrt (no domain checking) make sure to handle errors in calling method
_sqrt{T<:Float}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T, x)))

# x is double, z is double
@inline function sqrt(x::Double)
    r = _sqrt(x.hi)
    Double(r, (x.lo + fma(-r, r, x.hi))/(r+r))
end

# x is single, z is double
@inline function sqrt(x::Single)
    r = _sqrt(x.hi)
    Double(r, fma(-r, r, x.hi)/(r+r))
end

### auxiliary

==(x::Double, y::Single) = x.hi == y.hi
==(x::Single, y::Double) = y == x

abs(x::Single) = Single(abs(x.hi))
abs(x::Double) = Double(abs(x.hi), abs(x.lo))

scale{T<:Float}(x::Double{T}, s::T) = Double(s*x.hi, s*x.lo)
scale{T<:Float}(s::T, x::Double{T}) = Double(s*x.hi, s*x.lo)
scale{T<:Float}(x::Single{T}, s::T) = Single(s*x.hi)
scale{T<:Float}(s::T, x::Single{T}) = Single(s*x.hi)

function show{T}(io::IO, x::Double{T})
    println(io, "Double{$T}")
    print(io, x.hi, ", ", x.lo)
end

# hack :P
function show{T}(io::IO, x::Single{T})
    println(io, "Double{$T}")
    print(io, x.hi)
end


# Determine if hardware FMA is available, should probably check with LLVM, see #9855.
# Checks if the `fma` function is fast for the floating point type `T`: typically is it a
# native instruction (`true`) or does it fall back on a software implementation (`false`).
function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)),nextfloat(one(T)),-nextfloat(one(T),2)) != zero(T))
end
is_fma_fast() = is_fma_fast(Float64) && is_fma_fast(Float32)

if is_fma_fast()
    const twoprod = _twoprod_fma
    const pdiv = _pdiv_fma
else
    const twoprod = _twoprod
    const pdiv = _pdiv
end


# end