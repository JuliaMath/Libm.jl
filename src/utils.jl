"""
    highword(x::Union{Float32,Float64}) -> UInt32

Get the most significant 32 bits as a `UInt32` from `x`.
Corresponds to `GET_HIGH_WORD` in musl
"""
@inline highword(x::Float64) = (reinterpret(UInt64, x) >> 32) % UInt32
@inline highword(x::Float32) =  reinterpret(UInt32, x)


"""
    trunclo(x::Union{Float32,Float64})

Truncates the lower-order bits of `x` so that the result of the multiplication `trunclo(x)
* trunclo(y)` is exact, assuming no underflow or overflow occurs.
"""
@inline trunclo(x::Float64) =
    reinterpret(Float64, reinterpret(UInt64,x) & 0xffff_ffff_f800_0000) # truncate 27 bits
@inline trunclo(x::Float32) =
    reinterpret(Float32, reinterpret(UInt32,x) & 0xffff_f000)



"""
    setlowword(d::Float64, lo::UInt32)

Returns the least significant 32 bits of `d` to `lo`.
Corresponds to `SET_LOW_WORD` in musl
"""
@inline setlowword(d::Float64, lo::UInt32) = reinterpret(Float64, reinterpret(UInt64, d) & 0xffff_ffff_0000_0000 | lo)


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
