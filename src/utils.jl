"""
Get the more significant 32 bit int from a double.
Corresponds to GET_HIGH_WORD in musl
"""
highword(d::Float64) = (reinterpret(UInt64, d) >> 32) % UInt32

"""
Set the less significant 32 bits of the double d to the unsigned 32 bit integer lo
Corresponds to SET_LOW_WORD in musl
"""
lowword!(d::Float64, lo::UInt32) = reinterpret(Float64, reinterpret(UInt64, d) & 0xffffffff00000000 | lo)


# determine if hardware FMA is available
# should probably check with LLVM, see #9855.
"""
    is_fma_fast(T)

Checks if the `fma` function is fast for the floating point type `T`: typically is it a native instruction (`true`) or does it fall back on a software implementation (`false`).
    """
function is_fma_fast end
for T in (Float32, Float64)
    @eval is_fma_fast(::Type{$T}) = $(muladd(nextfloat(one(T)),nextfloat(one(T)),-nextfloat(one(T),2)) != zero(T))
end
