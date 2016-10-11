"""
    exp(x)

Compute the base ``e`` exponential of ``x``, in other words ``e^x``.
"""
function exp end

"""
    exp2(x)

Compute the base ``2`` exponential of ``x``, in other words ``2^x``.
"""
function exp2 end


#  Method
#    1. Argument reduction: Reduce x to an r so that |r| <= 0.5*ln2. Given x,
#       find r and integer k such that
#                x = k*log2(e) + r,  |r| <= 0.5*ln2.
#       Here r is represented as r = hi-lo for better accuracy.
#
#    2. Approximate exp(r) by a special rational function on [0, 0.5*ln2]:
#           R(r^2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r^4/360 + ...
#
#      A special Remez algorithm on [0, 0.5*ln2] is used to generate a
#       polynomial to approximate R. In other words,
#
#           R(z) ~ 2.0 + P1*z + P2*z^2 + P3*z^3 + P4*z^4 + P5*z^5,
#
#       where z=r*r.
# 
#       The computation of exp(r) thus becomes
#                               2*r
#               exp(r) = 1 + ----------
#                             R(r) - r
#                                  r*c(r)
#                      = 1 + r + ----------- (for better accuracy)
#                                 2 - c(r)
#       where
#               c(r) = r - (P1*r^2  + P2*r^4  + ... + P5*r^10 + ...).
#
#    3. Scale back: exp(x) = 2^k * exp(r)


MAXEXP(::Type{Float64}) = 7.09782712893383996732e2 # log 2^1023*(2-2^-52)
MAXEXP(::Type{Float32}) = 88.72283905206835f0      # log 2^127 *(2-2^-23)
MAXEXP(::Type{Float16}) = Float16(11.09)           # log 2^15  *(2-2^-10)

MAXEXP2(::Type{Float64}) = 1024         # log2 2^1023*(2-2^-52)
MAXEXP2(::Type{Float32}) = 128          # log2 2^127 *(2-2^-23)
MAXEXP2(::Type{Float16}) = Float16(16)  # log2 2^15  *(2-2^-10)

# one less than the min exponent since we can sqeeze a bit more from the exp function
MINEXP(::Type{Float64}) = -7.451332191019412076235e2                     # log 2^-1075
MINEXP(::Type{Float32}) = -103.972077083991796412584818218726485211325f0 # log 2^-150
MINEXP(::Type{Float16}) = Float16(-17.32868)                             # log 2^-25

MINEXP2(::Type{Float64}) = -1075        # log 2^-1075
MINEXP2(::Type{Float32}) = -150         # log 2^-150
MINEXP2(::Type{Float16}) = Float16(-25) # log 2^-25


# coefficients from:
# origin: FreeBSD /usr/src/lib/msun/src/e_exp.c */
# ====================================================
# Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
# Permission to use, copy, modify, and distribute this
# software is freely granted, provided that this notice
# is preserved.
# ====================================================
@inline @oftype function _exp{T}(hi::T, lo::T)
    r = hi - lo
    z = r*r
    p = r - z *
    (1.66666666666666019037e-01  + z *
    (-2.77777777770155933842e-03 + z *
    (6.61375632143793436117e-05  + z *
    (-1.65339022054652515390e-06 + z *
    4.13813679705723846039e-08))))
    return 1.0 - ((lo - (r*p)/(2.0 - p)) - hi)
end

# custom coefficients
@inline @oftype function _exp{T<:SmallFloatTypes}(hi::T, lo::T)
    r = hi - lo
    z = r*r
    p = r - z *
    (0.1666666567325592041015625 + z *
    (-2.777527086436748504638671875e-3 + z *
    6.451140507124364376068115234375e-5))
    return 1.0 - ((lo - (r*p)/(2.0 - p)) - hi)
end


function exp{T}(x::T)
    x > MAXEXP(T) && return T(Inf)
    x < MINEXP(T) && return T(0.0)
 
    # reduce; computed as r = hi - lo for extra precision.
    k = round(T(LOG2E)*x) # k is a float type
    hi = muladd(k, -LN2U(T), x)
    lo = k*LN2L(T)

    # compute approximation
    x = _exp(hi,lo)
    n = _trunc(k) # ldexp needs n as an int
    return _ldexp(x,n)
end


#  Method: we simply scale the input argument by log(2) and use the
#  coefficients from exp(x)
function exp2{T}(x::T)
    x > MAXEXP2(T) && return T(Inf)
    x < MINEXP2(T) && return T(0)
 
    # reduce; computed as r = hi - lo for extra precision.
    k = round(x)
    t = x - k
    hi = t*LN2U(T)
    lo = -t*LN2L(T)

    # compute approximation
    x = _exp(hi,lo)
    n = _trunc(k) # ldexp needs n as an int
    return _ldexp(x,n)
end
