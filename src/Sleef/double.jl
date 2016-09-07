immutable Double{T<:FloatTypes}
    hi::T
    lo::T
end
Double{T}(x::T) = Double(x, zero(T))

@inline trunclo(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xffff_ffff_f800_0000) # clear lower 27 bits (leave upper 26 bits)
@inline trunclo(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xffff_f000) # clear lowest 12 bits (leave upper 12 bits)

@inline function splitprec(x::FloatTypes)
    hx = trunclo(x)
    hx, x-hx
end

@inline function ddnormalize{T}(x::Double{T})
    r = x.hi + x.lo
    Double(r, (x.hi - r) + x.lo)
end

@inline ddscale{T<:FloatTypes}(x::Double{T}, s::T) = Double(s*x.hi, s*x.lo)

@inline ddneg{T}(x::Double{T}) = Double(-x.hi,-x.lo)

# quick-two-sum x+y
@inline function ddadd{T<:FloatTypes}(x::T, y::T) #WARNING |x| >= |y|
    r = x + y
    Double(r, (x - r) + y)
end

@inline function ddadd{T<:FloatTypes}(x::T, y::Double{T}) #WARNING |x| >= |y|
    r = x + y.hi
    Double(r, (x - r) + y.hi + y.lo)
end

@inline function ddadd{T<:FloatTypes}(x::Double{T}, y::T) #WARNING |x| >= |y|
    r = x.hi + y
    Double(r, (x.hi - r) + y + x.lo)
end

@inline function ddadd{T}(x::Double{T}, y::Double{T}) #WARNING |x| >= |y|
    r = x.hi + y.hi
    Double(r, (x.hi - r) + y.hi + y.lo + x.lo)
end

@inline function ddsub{T}(x::Double{T}, y::Double{T}) #WARNING |x| >= |y|
    r = x.hi - y.hi
    Double(r, (x.hi - r) - y.hi - y.lo + x.lo)
end

# two-sum x+y  NO BRANCH
@inline function ddadd2{T<:FloatTypes}(x::T, y::T)
    s = x + y
    v = s - x
    Double(s, (x - (s - v)) + (y - v))
end

@inline function ddadd2{T<:FloatTypes}(x::T, y::Double{T})
    s = x + y.hi
    v = s - x
    Double(s, (x - (s - v)) + (y.hi - v) + y.lo)
end

@inline ddadd2{T<:FloatTypes}(x::Double{T}, y::T) = ddadd2(y,x)

@inline function ddadd2{T}(x::Double{T}, y::Double{T})
    s = x.hi + y.hi
    v = s - x.hi
    Double(s,(x.hi - (s - v)) + (y.hi - v) + x.lo + y.lo)
end

if is_fma_fast()

    # two-prod-fma
    @inline function ddmul{T<:FloatTypes}(x::T, y::T)
        z = x*y
        Double(z, fma(x, y, -z))
    end

    @inline function ddmul{T<:FloatTypes}(x::Double{T}, y::T)
        z = x.hi*y
        Double(z, fma(x.hi, y, -z) + x.lo*y)
    end

    @inline ddmul{T<:FloatTypes}(x::T, y::Double{T}) = ddmul(y,x)

    @inline function ddmul{T}(x::Double{T}, y::Double{T})
        z = x.hi*y.hi
        Double(z, fma(x.hi, y.hi, -z) + x.hi*y.lo + x.lo*y.hi)
    end

    # x^2
    @inline function ddsqu{T<:FloatTypes}(x::T)
        z = x*x
        Double(z, fma(x,x,-z))
    end

    @inline function ddsqu{T}(x::Double{T})
        z = x.hi*x.hi
        Double(z, fma(x.hi, x.hi, -z) + x.hi*(x.lo+x.lo))
    end

else

    #two-prod x*y
    @inline function ddmul{T<:FloatTypes}(x::T, y::T)
        hx, lx = splitprec(x)
        hy, ly = splitprec(y)
        z = x*y
        Double(z, ((hx*hy-z) + lx*hy + hx*ly) + lx*ly)
    end

    @inline function ddmul{T<:FloatTypes}(x::Double{T}, y::T)
        hx, lx = splitprec(x.hi)
        hy, ly = splitprec(y)
        z = x.hi*y
        Double(z, (hx*hy-z) + lx*hy + hx*ly + lx*ly + x.lo*y)
    end

    @inline ddmul{T<:FloatTypes}(x::T, y::Double{T}) = ddmul(y,x)

    @inline function ddmul{T}(x::Double{T}, y::Double{T})
        hx, lx = splitprec(x.hi)
        hy, ly = splitprec(y.hi)
        z = x.hi*y.hi
        Double(z, (((hx*hy-z) + lx*hy + hx*ly) + lx*ly) + x.hi*y.lo + x.lo*y.hi)
    end

    # x^2
    @inline function ddsqu{T<:FloatTypes}(x::T)
        hx, lx = splitprec(x)
        z = x*x
        Double(z, (hx*hx-z) + lx*(hx + hx) + lx*lx)
    end

    @inline function ddsqu{T}(x::Double{T})
        hx, lx = splitprec(x.hi)
        z = x.hi*x.hi
        Double(z, (hx*hx-z) + lx*(hx+hx) + lx*lx + x.hi*(x.lo+x.lo))
    end

end

# x/y
@inline function dddiv{T}(x::Double{T}, y::Double{T})
    invy = 1/y.hi
    c = x.hi*invy
    u = ddmul(c, y.hi)
    Double(c,((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)*invy)
end

# 1/x
@inline function ddrec{T<:FloatTypes}(x::T)
    invx = 1/x
    c = one(T)*invx
    u = ddmul(c,x)
    Double(c, (one(T) - u.hi - u.lo)*invx)
end

@inline function ddrec{T}(x::Double{T})
    invx = 1/x.hi
    c = one(T)*invx
    u = ddmul(c,x.hi)
    Double(c, (one(T) -u.hi - u.lo - c*x.lo)*invx)
end

# sqrt(x)
@inline function ddsqrt{T}(x::Double{T})
    c = _sqrt(x.hi)
    u = ddsqu(c)
    Double(c, (x.hi - u.hi - u.lo + x.lo)*T(0.5)/c)
end
