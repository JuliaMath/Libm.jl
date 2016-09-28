import Base: -, flipsign, normalize, scale, convert, <

immutable Double{T<:FloatTypes} <: Number
    hi::T
    lo::T
end
Double{T}(x::T) = Double(x, zero(T))

convert{T<:FloatTypes}(::Type{T}, x::Double) = x.hi + x.lo


@inline trunclo(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xffff_ffff_f800_0000) # clear lower 27 bits (leave upper 26 bits)
@inline trunclo(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xffff_f000) # clear lowest 12 bits (leave upper 12 bits)

@inline function splitprec(x::FloatTypes)
    hx = trunclo(x)
    hx, x-hx
end

@inline function normalize{T}(x::Double{T})
    r = x.hi + x.lo
    Double(r, (x.hi - r) + x.lo)
end

@inline flipsign{T<:FloatTypes}(x::Double{T}, y::T) = Double(flipsign(x.hi, y), flipsign(x.lo, y))

@inline scale{T<:FloatTypes}(x::Double{T}, s::T) = Double(s*x.hi, s*x.lo)

@inline -{T}(x::Double{T}) = Double(-x.hi,-x.lo)

@inline function <{T}(x::Double{T}, y::Double{T})
    x.hi < y.hi
end

@inline function <{T}(x::Double{T}, y::Number)
    x.hi < y
end

@inline function <{T}(x::Number, y::Double{T})
    x < y.hi
end

# quick-two-sum x+y
@inline function dadd{T<:FloatTypes}(x::T, y::T) #WARNING |x| >= |y|
    s = x + y
    Double(s, (x - s) + y)
end

@inline function dadd{T<:FloatTypes}(x::T, y::Double{T}) #WARNING |x| >= |y|
    s = x + y.hi
    Double(s, (x - s) + y.hi + y.lo)
end

@inline function dadd{T<:FloatTypes}(x::Double{T}, y::T) #WARNING |x| >= |y|
    s = x.hi + y
    Double(s, (x.hi - s) + y + x.lo)
end

@inline function dadd{T}(x::Double{T}, y::Double{T}) #WARNING |x| >= |y|
    s = x.hi + y.hi
    Double(s, (x.hi - s) + y.hi + y.lo + x.lo)
end

@inline function dsub{T}(x::Double{T}, y::Double{T}) #WARNING |x| >= |y|
    s = x.hi - y.hi
    Double(s, (x.hi - s) - y.hi - y.lo + x.lo)
end

@inline function dsub{T}(x::Double{T}, y::T) #WARNING |x| >= |y|
    s = x.hi - y
    Double(s, (x.hi - s) - y + x.lo)
end

@inline function dsub{T}(x::T, y::Double{T}) #WARNING |x| >= |y|
    s = x - y.hi
    Double(s, (x - s) - y.hi - y.lo)
end

@inline function dsub{T}(x::T, y::T) #WARNING |x| >= |y|
    s = x - y
    Double(s, (x - s) - y)
end


# two-sum x+y  NO BRANCH
@inline function dadd2{T<:FloatTypes}(x::T, y::T)
    s = x + y
    v = s - x
    Double(s, (x - (s - v)) + (y - v))
end

@inline function dadd2{T<:FloatTypes}(x::T, y::Double{T})
    s = x + y.hi
    v = s - x
    Double(s, (x - (s - v)) + (y.hi - v) + y.lo)
end

@inline dadd2{T<:FloatTypes}(x::Double{T}, y::T) = dadd2(y,x)

@inline function dadd2{T}(x::Double{T}, y::Double{T})
    s = x.hi + y.hi
    v = s - x.hi
    Double(s,(x.hi - (s - v)) + (y.hi - v) + x.lo + y.lo)
end

@inline function dsub2{T<:FloatTypes}(x::T, y::T)
    s = x - y
    v = s - x
    Double(s, (x - (s - v)) + (-y - v))
end

@inline function dsub2{T<:FloatTypes}(x::T, y::Double{T})
    s = x - y.hi
    v = s - x
    Double(s, (x - (s - v)) + (-y.hi - v) - y.lo)
end

@inline dsub2{T<:FloatTypes}(x::Double{T}, y::T) = dsub2(y,-x)

@inline function dsub2{T}(x::Double{T}, y::Double{T})
    s = x.hi - y.hi
    v = s - x.hi
    Double(s,(x.hi - (s - v)) + (-y.hi - v) + x.lo - y.lo)
end



if is_fma_fast()

    # two-prod-fma
    @inline function dmul{T<:FloatTypes}(x::T, y::T)
        z = x*y
        Double(z, fma(x, y, -z))
    end

    @inline function dmul{T<:FloatTypes}(x::Double{T}, y::T)
        z = x.hi*y
        Double(z, fma(x.hi, y, -z) + x.lo*y)
    end

    @inline dmul{T<:FloatTypes}(x::T, y::Double{T}) = dmul(y,x)

    @inline function dmul{T}(x::Double{T}, y::Double{T})
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

    # sqrt(x)
    @inline function dsqrt{T}(x::Double{T})
        zhi = _sqrt(x.hi)
        Double(zhi, (x.lo + fma(-zhi, zhi, x.hi))/(zhi+zhi))
    end

    # x/y
    @inline function ddiv{T}(x::Double{T}, y::Double{T})
        invy = 1/y.hi
        zhi = x.hi*invy
        Double(zhi, (fma(-zhi, y.hi, x.hi) + fma(-zhi, y.lo, x.lo))*invy)
    end

    # 1/x
    @inline function ddrec{T<:FloatTypes}(x::T)
        zhi = 1/x
        Double(zhi, fma(-zhi, x, one(T))*zhi)
    end

    @inline function ddrec{T}(x::Double{T})
        zhi = 1/x.hi
        Double(zhi, (fma(-zhi, x.hi, one(T)) + -zhi*x.lo)*zhi)
    end

else

    #two-prod x*y
    @inline function dmul{T<:FloatTypes}(x::T, y::T)
        hx, lx = splitprec(x)
        hy, ly = splitprec(y)
        z = x*y
        Double(z, ((hx*hy-z) + lx*hy + hx*ly) + lx*ly)
    end

    @inline function dmul{T<:FloatTypes}(x::Double{T}, y::T)
        hx, lx = splitprec(x.hi)
        hy, ly = splitprec(y)
        z = x.hi*y
        Double(z, (hx*hy-z) + lx*hy + hx*ly + lx*ly + x.lo*y)
    end

    @inline dmul{T<:FloatTypes}(x::T, y::Double{T}) = dmul(y,x)

    @inline function dmul{T}(x::Double{T}, y::Double{T})
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

    # sqrt(x)
    @inline function dsqrt{T}(x::Double{T})
        c = _sqrt(x.hi)
        u = ddsqu(c)
        Double(c, (x.hi - u.hi - u.lo + x.lo)/(c+c))
    end

    # x/y
    @inline function ddiv{T}(x::Double{T}, y::Double{T})
        invy = 1/y.hi
        c = x.hi*invy
        u = dmul(c, y.hi)
        Double(c,((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)*invy)
    end

    # 1/x
    @inline function ddrec{T<:FloatTypes}(x::T)
        c = 1/x
        u = dmul(c,x)
        Double(c, (one(T) - u.hi - u.lo)*c)
    end

    @inline function ddrec{T}(x::Double{T})
        c = 1/x.hi
        u = dmul(c,x.hi)
        Double(c, (one(T) - u.hi - u.lo - c*x.lo)*c)
    end

end




