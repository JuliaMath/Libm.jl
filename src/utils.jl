# Get the more significant 32 bit int from a double
function get_high_word(d::Float64)
    u = reinterpret(UInt64, d)
    return (u>>32) % UInt32
end

# Set the less significant 32 bits of a double from an int
function set_low_word(d::Float64, lo::UInt32)
    return reinterpret(Float64, reinterpret(UInt64, d) & 0xffffffff00000000 | lo)
end

# Scale number x*2^n
function scalbn(x::Float64, n::Integer)
    if n > 1023
        x *= 0x1p1023
        n -= 1023;
        if n > 1023
            x *= 0x1p1023
            n -= 1023;
            if n > 1023
                n = 1023
             end
        end
    elseif n < -1022
        x *= 0x1p-1022
        n += 1022
        if n < -1022
            x *= 0x1p-1022
            n += 1022
            if n < -1022
                n = -1022
            end
        end
    end
    u = reinterpret(Float64, (0x3ff+n) % UInt64 << 52)
    return x*u
end
