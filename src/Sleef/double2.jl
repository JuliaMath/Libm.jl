
immutable Double2{T<:FTypes}
    x::T
    y::T
end

# NDEBUG = true
# if NDEBUG
#     function checkfp(x::FTypes)
#       isinf(x) || isnan(x) && return true
#       return false
#     end
# end

upper(d::Float64) = reinterpret(Float64, reinterpret(UInt64, d) & 0xfffffffff8000000)
upper(d::Float32) = reinterpret(Float32, reinterpret(UInt32, d) & 0xfffff000)

function ddnormalize_d2_d2(t::Double2)
    sx = t.x + t.y
    sy = t.x - sx + t.y
    return Double2(sx,sy)
end

function ddscale_d2_d2_d{T<:FTypes}(d::Double2{T}, s::T)
    rx = d.x * s
    ry = d.y * s
    return Double2(rx,ry)
end

function ddneg_d2_d2(d::Double2)
    return Double2(-d.x,-d.y)
end

function ddadd_d2_d_d{T<:FTypes}(x::T, y::T)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x) || checkfp(y) || abs(x) >= abs(y)) && error(@sprintf "[ddadd_d2_d_d : %g, %g]" x y)
    # end
    rx = x + y
    ry = x - rx + y
  return Double2(rx,ry)
end

function ddadd2_d2_d_d{T<:FTypes}(x::T, y::T)
    rx = x + y
    v = rx - x
    ry = (x - (rx - v)) + (y - v)
    return Double2(rx,ry)
end

function ddadd_d2_d2_d{T<:FTypes}(x::Double2{T}, y::T)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y) || abs(x.x) >= abs(y)) && error(@sprintf "[ddadd_d2_d2_d : %g %g]" x.x y)
    # end
    rx = x.x + y
    ry = x.x - rx + y + x.y
    return Double2(rx,ry)
end

function ddadd2_d2_d2_d{T<:FTypes}(x::Double2{T}, y::T)
    # |x| >= |y|
    rx  = x.x + y;
    v = rx - x.x;
    ry = (x.x - (rx - v)) + (y - v);
    ry += x.y;
    return Double2(rx,ry)
end

function ddadd_d2_d_d2{T<:FTypes}(x::T, y::Double2{T})
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x) || checkfp(y.x) || abs(x) >= abs(y.x)) && error(@sprintf "[ddadd_d2_d_d2 : %g %g]" x y.x)
    # end
    rx = x + y.x;
    ry = x - rx + y.x + y.y
    return Double2(rx,ry)
end

function ddadd2_d2_d_d2{T<:FTypes}(x::T, y::Double2{T})
    rx  = x + y.x
    v = rx - x
    ry = (x - (rx - v)) + (y.x - v) + y.y
    return Double2(rx,ry)
end

function ddadd_d2_d2_d2(x::Double2, y::Double2)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y.x) || abs(x.x) >= abs(y.x)) && error(@sprintf "[ddadd_d2_d2_d2 : %g %g]" x.x y.x)
    # end
    rx = x.x + y.x
    ry = x.x - rx + y.x + x.y + y.y
  return Double2(rx,ry)
end

function ddadd2_d2_d2_d2(x::Double2, y::Double2)
    rx  = x.x + y.x
    v = rx - x.x
    ry = (x.x - (rx - v)) + (y.x - v)
    ry += x.y + y.y
    return Double2(rx,ry)
end

function ddsub_d2_d2_d2(x::Double2, y::Double2)
    # |x| >= |y|
    # if NDEBUG
    #     !(checkfp(x.x) || checkfp(y.x) || abs(x.x) >= abs(y.x)) && error(@sprintf "[ddsub_d2_d2_d2 : %g %g]" x.x y.x)
    # end
    rx = x.x - y.x;
    ry = x.x - rx - y.x + x.y - y.y;
    return Double2(rx,ry)
end

function dddiv_d2_d2_d2{T<:FTypes}(n::Double2{T}, d::Double2{T})
    t = one(T)/d.x
    dh  = upper(d.x); dl  = d.x - dh
    th  = upper(t);   tl  = t   - th
    nhh = upper(n.x); nhl = n.x - nhh
    qx = n.x * t;
    u = -qx + nhh * th + nhh * tl + nhl * th + nhl * tl +
            qx * (one(T) - dh * th - dh * tl - dl * th - dl * tl)
    qy = t * (n.y - qx * d.y) + u
    return Double2(qx,qy)
end


function ddmul_d2_d_d{T<:FTypes}(x::T, y::T)
    xh = upper(x); xl = x - xh
    yh = upper(y); yl = y - yh
    rx = x * y
    ry = xh * yh - rx + xl * yh + xh * yl + xl * yl
    return Double2(rx,ry)
end

function ddmul_d2_d2_d{T<:FTypes}(x::Double2{T}, y::T)
    xh = upper(x.x); xl = x.x - xh
    yh = upper(y);   yl = y   - yh
    rx = x.x * y
    ry = xh * yh - rx + xl * yh + xh * yl + xl * yl + x.y * y
    return Double2(rx,ry)
end

function ddmul_d2_d2_d2(x::Double2, y::Double2)
    xh = upper(x.x); xl = x.x - xh
    yh = upper(y.x); yl = y.x - yh
    rx = x.x * y.x
    ry = xh * yh - rx + xl * yh + xh * yl + xl * yl + x.x * y.y + x.y * y.x
    return Double2(rx,ry)
end

function ddsqu_d2_d2(x::Double2)
    xh = upper(x.x); xl = x.x - xh
    rx = x.x * x.x
    ry = xh * xh - rx + (xh + xh) * xl + xl * xl + x.x * (x.y + x.y)
    return Double2(rx,ry)
end

function ddrec_d2_d{T<:FTypes}(d::T)
    t = one(T)/d
    dh = upper(d); dl = d - dh
    th = upper(t); tl = t - th
    qx = t
    qy = t * (one(T) - dh * th - dh * tl - dl * th - dl * tl)
    return Double2(qx,qy)
end

function ddrec_d2_d2{T<:FTypes}(d::Double2{T})
    t = one(T)/d.x
    dh = upper(d.x); dl = d.x - dh
    th = upper(t);   tl = t   - th
    qx = t
    qy = t * (one(T) - dh * th - dh * tl - dl * th - dl * tl - d.y * t)
    return Double2(qx,qy)
end

function ddsqrt_d2_d2{T<:FTypes}(d::Double2{T})
    t = sqrt(d.x + d.y)
    return ddscale_d2_d2_d(ddmul_d2_d2_d2(ddadd2_d2_d2_d2(d, ddmul_d2_d_d(t, t)), ddrec_d2_d(t)), T(0.5))
end
