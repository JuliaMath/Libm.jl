#Relates to OpenLibm: src/math_private, however little of that code can be transfered as it uses C unions.
using Base: @pure

"""Return the least significant word of d.
Coresponds to OpenLibm GET_LOW_WORD(ret,d)"""
@inline @pure lowword(d::Float64) = reinterpret(UInt64, d) % UInt32

"""Returns the most significant word of d
Coresponds to OpenLibm GET_HIGH_WORD(ret,d)
"""
@inline @pure highword(d::Float64) = (reinterpret(UInt64, d) >> 32) % UInt32 

"""
Combined the high and low words, to make a Float64.
Coresponds to OpenLibm INSERTWORDS(ret, hw, lw)
"""
@inline @pure function combinewords(hw::UInt32, lw::UInt32)::Float64
	reinterpret(Float64, (hw % UInt64) << 32 + lw)
end
