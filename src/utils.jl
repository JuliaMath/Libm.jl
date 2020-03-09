
# Useful utilities from Base
using Base.Math: sign_mask, exponent_mask, significand_mask, exponent_one, exponent_bias, significand_bits, @horner

# Similar to @horner, but converts coefficients to same type as x
# TODO this seems to be broken for certain type conversions
macro horner_oftype(x, p...)
    ex = :(oftype($x,$(esc(p[end]))))
    for i = length(p)-1:-1:1
        ex = :(muladd(t, $ex, oftype($x,$(esc(p[i])))))
    end
    Expr(:block, :(t = $(esc(x))), ex)
end


"""
    highword(x::Union{Float32,Float64}) -> UInt32

Get the most significant 32 bits as a `UInt32` from `x`.
Corresponds to `GET_HIGH_WORD` in musl
"""
@inline highword(x::Float64) = UInt32(reinterpret(UInt64, x) >> 32)
@inline highword(x::Float32) =  reinterpret(UInt32, x)


"""
    trunclo(x::Union{Float32,Float64})

Truncates the lower order bits of `x` so that the result of the multiplication `trunclo(x)
* trunclo(y)` is exact, assuming no underflow or overflow occurs.

This relies on the following property: if `a` has `n` significant bits, and `b` has `m`
significant bits, then exact product `a*b` has either `n+m-1` or `n+m` significant bits:

* `Float64` has 53 significant bits (including implicit leading bit), so we need to
  truncate the last 27 significant bits (leaving 26 bits).

* `Float32` has 24 significant bits (including implicit leading bit), so we need to
  truncate the last 12 significant bits (leaving 12 bits).

This is typically faster than other methods of truncating lower order bits (such as
Veltkamp splitting, or converting `Float64`s to `Float32`s and back again). For LLVM 3.8
or greater, this should give optimal `ANDPD`/`ANDSD` instructions on supported x86
architectures, which doesn't require moving registers (Julia issue #9868).

NOTE: For odd significand sizes, such as `Float64`, when used as a replacement Veltkamp
splitting for computing extended precision multiplications with a Dekker-style `mul12`
algorithm, it can lose the last bit of precision. For example the function:

    function incorrect_mul2(x,y)
        hx = trunclo(x); tx = x-hx
        hy = trunclo(y); ty = y-hy
        p = hx*hy
        q = hx*ty + tx*hy
        z = p+q
        zz = p-z+q+tx*ty
        z, zz
    end

will give an incorrect result for the case `x = y = 0x1.800000e000001p+0`.
"""
@inline trunclo(x::Float64) =
    reinterpret(Float64, reinterpret(UInt64,x) & 0xffff_ffff_f800_0000)
@inline trunclo(x::Float32) =
    reinterpret(Float32, reinterpret(UInt32,x) & 0xffff_f000)


# determine if hardware FMA is available
# should probably check with LLVM, see #9855.
"""
    is_fma_fast(T)

Checks if the `fma` function is fast for the floating point type `T`: typically is it a
native instruction (`true`) or does it fall back on a software implementation (`false`).
"""
function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)),nextfloat(one(T)),-nextfloat(one(T),2)) != zero(T))
end
