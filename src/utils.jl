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

