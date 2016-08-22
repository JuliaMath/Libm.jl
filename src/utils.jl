"""
Get the more significant 32 bit int from a double.
Corresponds to GET_HIGH_WORD in musl
"""
highword(d::Float64) = (reinterpret(UInt64, d) >> 32) % UInt32

# Scale number x*2^n
function scalbn(x::Float64, n::Int32)
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