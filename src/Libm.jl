module Libm

include("erf.jl")


# Utils

# Get the more significant 32 bit int from a double
function GET_HIGH_WORD(d::Float64)
    u = reinterpret(UInt64, d)
    (u >> 32) % UInt32
end

# Set the less significant 32 bits of a double from an int
SET_LOW_WORD(d::Float64, lo::UInt32) = reinterpret(Float64, reinterpret(UInt64, d) & 0xffffffff00000000 | lo)

end
