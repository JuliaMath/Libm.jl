# Get the more significant 32 bit int from a double
function get_high_word(d::Float64)
    u = reinterpret(UInt64, d)
    (u>>32)%UInt32
end

# Set the less significant 32 bits of a double from an int
set_low_word(d::Float64, lo::UInt32) = reinterpret(Float64, reinterpret(UInt64, d) & 0xffffffff00000000 | lo)

