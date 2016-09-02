
immutable Double{T<:FloatTypes}
    x::T
    y::T
end

# NDEBUG = true
# if NDEBUG
#     function checkfp(x::FloatTypes)
#       isinf(x) || isnan(x) && return true
#       return false
#     end
# end

@inline function ddnormalize_d2_d2(t::Double)
    sx = t.x + t.y
    sy = t.x - sx + t.y
    return Double(sx,sy)
end

@inline function ddscale_d2_d2_d{T<:FloatTypes}(d::Double{T}, s::T)
    return Double(d.x * s, d.y * s)
end

@inline function ddneg_d2_d2(d::Double)
    return Double(-d.x,-d.y)
end

# quick-two-sum
@inline function ddadd_d2_d_d{T<:FloatTypes}(x::T, y::T)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x) || checkfp(y) || abs(x) >= abs(y)) && error(@sprintf "[ddadd_d2_d_d : %g, %g]" x y)
    # end
    rx = x + y
    ry = y - (rx - x)
  return Double(rx,ry)
end

# two-sum
@inline function ddadd2_d2_d_d{T<:FloatTypes}(x::T, y::T)
    rx = x + y
    v = rx - x
    ry = (x - (rx - v)) + (y - v)
    return Double(rx,ry)
end

@inline function ddadd_d2_d2_d{T<:FloatTypes}(x::Double{T}, y::T)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y) || abs(x.x) >= abs(y)) && error(@sprintf "[ddadd_d2_d2_d : %g %g]" x.x y)
    # end
    rx = x.x + y
    ry = x.x - rx + y + x.y
    return Double(rx,ry)
end

@inline function ddadd2_d2_d2_d{T<:FloatTypes}(x::Double{T}, y::T)
    # |x| >= |y|
    rx  = x.x + y;
    v = rx - x.x;
    ry = (x.x - (rx - v)) + (y - v);
    ry += x.y;
    return Double(rx,ry)
end

@inline function ddadd_d2_d_d2{T<:FloatTypes}(x::T, y::Double{T})
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x) || checkfp(y.x) || abs(x) >= abs(y.x)) && error(@sprintf "[ddadd_d2_d_d2 : %g %g]" x y.x)
    # end
    rx = x + y.x;
    ry = x - rx + y.x + y.y
    return Double(rx,ry)
end

@inline function ddadd2_d2_d_d2{T<:FloatTypes}(x::T, y::Double{T})
    rx  = x + y.x
    v = rx - x
    ry = (x - (rx - v)) + (y.x - v) + y.y
    return Double(rx,ry)
end

@inline function ddadd_d2_d2_d2(x::Double, y::Double)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y.x) || abs(x.x) >= abs(y.x)) && error(@sprintf "[ddadd_d2_d2_d2 : %g %g]" x.x y.x)
    # end
    rx = x.x + y.x
    ry = x.x - rx + y.x + x.y + y.y
  return Double(rx,ry)
end

@inline function ddadd2_d2_d2_d2(x::Double, y::Double)
    rx  = x.x + y.x
    v = rx - x.x
    ry = (x.x - (rx - v)) + (y.x - v)
    ry += x.y + y.y
    return Double(rx,ry)
end

@inline function ddsub_d2_d2_d2(x::Double, y::Double)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y.x) || abs(x.x) >= abs(y.x)) && error(@sprintf "[ddsub_d2_d2_d2 : %g %g]" x.x y.x)
    # end
    rx = x.x - y.x;
    ry = x.x - rx - y.x + x.y - y.y;
    return Double(rx,ry)
end

@inline function dddiv_d2_d2_d2{T<:FloatTypes}(n::Double{T}, d::Double{T})
    t = one(T)/d.x
    dh  = upper(d.x); dl  = d.x - dh
    th  = upper(t);   tl  = t   - th
    nhh = upper(n.x); nhl = n.x - nhh
    qx = n.x * t;
    u = -qx + nhh * th + nhh * tl + nhl * th + nhl * tl +
            qx * (one(T) - dh * th - dh * tl - dl * th - dl * tl)
    qy = t * (n.y - qx * d.y) + u
    return Double(qx,qy)
end

# two-prod
@inline function ddmul_d2_d_d{T<:FloatTypes}(x::T, y::T)
    hx = upper(x); lx = x - hx
    hy = upper(y); ly = y - hy
    rx = x*y
    ry = ((hx*hy - rx) + lx*hy + hx*ly) + lx*ly
    return Double(rx,ry)
end

# two-prod-fma
@inline function ddmul_d2_d_d_fma{T<:FloatTypes}(x::T, y::T)
    rx = x*y
    ry = muladd(x, y, -rx)
    return Double(rx,ry)
end

@inline function ddmul_d2_d2_d{T<:FloatTypes}(x::Double{T}, y::T)
    hx = upper(x.x); lx = x.x - hx
    hy = upper(y);   ly = y   - hy
    rx = x.x * y
    ry = hx * hy - rx + lx * hy + hx * ly + lx * ly + x.y * y
    return Double(rx,ry)
end

@inline function ddmul_d2_d2_d2(x::Double, y::Double)
    hx = upper(x.x); lx = x.x - hx
    hy = upper(y.x); ly = y.x - hy
    rx = x.x * y.x
    ry = hx * hy - rx + lx * hy + hx * ly + lx * ly + x.x * y.y + x.y * y.x
    return Double(rx,ry)
end

@inline function ddsqu_d2_d2(x::Double)
    hx = upper(x.x); lx = x.x - hx
    rx = x.x * x.x
    ry = hx * hx - rx + (hx + hx) * lx + lx * lx + x.x * (x.y + x.y)
    return Double(rx,ry)
end

@inline function ddrec_d2_d{T<:FloatTypes}(d::T)
    t = one(T)/d
    dh = upper(d); dl = d - dh
    th = upper(t); tl = t - th
    qx = t
    qy = t * (one(T) - dh * th - dh * tl - dl * th - dl * tl)
    return Double(qx,qy)
end

@inline function ddrec_d2_d2{T<:FloatTypes}(d::Double{T})
    t = one(T)/d.x
    dh = upper(d.x); dl = d.x - dh
    th = upper(t);   tl = t   - th
    qx = t
    qy = t * (one(T) - dh * th - dh * tl - dl * th - dl * tl - d.y * t)
    return Double(qx,qy)
end

@inline function ddsqrt_d2_d2{T<:FloatTypes}(d::Double{T})
    t = _sqrt(d.x + d.y)
    return ddscale_d2_d2_d(ddmul_d2_d2_d2(ddadd2_d2_d2_d2(d, ddmul_d2_d_d(t, t)), ddrec_d2_d(t)), T(0.5))
end
