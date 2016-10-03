# Float128
# C1 + C2 = ln 2 
C1{T<:Float}(::Type{T}) = T(6.93145751953125e-1)
C2{T<:Float}(::Type{T}) = T(1.428606820309417232121458176568075500134e-6)

C1{T<:SFloat}(::Type{T}) =  T(0.693359375)
C2{T<:SFloat}(::Type{T}) =  T(-2.12194440e-4)

let

# FreeBSD's coefficients /usr/src/lib/msun/src/e_exp.c
#   Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
#   Permission to use, copy, modify, and distribute this software is freely
#   granted, provided that this notice is preserved.
const a1 = 1.66666666666666019037e-01 
const a2 = -2.77777777770155933842e-03
const a3 = 6.61375632143793436117e-05 
const a4 = -1.65339022054652515390e-06
const a5 = 4.13813679705723846039e-08 

# custom coefficients @musm. TODO: just build a polynomial
const b1 = 1.666666499273473230749034e-01 
const b2 = -2.777534742424458760198523e-03
const b3 = 6.500773434631342911309326e-05 

global @inline function _exp{T}(xh::T,xl::T)
    x = xh - xl
    z = x*x  
    px = x - z*(@horner_oftype z a1 a2 a3 a4 a5)
    return T(1.0) - ((xl - (x*px)/(T(2.0)-px)) - xh)
end

global @inline function _exp{T<:SFloat}(xh::T,xl::T)
    x = xh - xl
    z = x*x  
    px = x - z*(@horner_oftype z b1 b2 b3)
    return T(1.0) - ((xl - (x*px)/(T(2.0)-px)) - xh)
end

global function exp{T}(x::T)
    isnan(x) && return x
    x > MAXLOG(T) && return T(Inf)
    x < MINLOG(T) && return T(0)
 
    px = round(T(LOG2E)*x)
    n = unsafe_trunc(Int,px)
    
    xh = muladd(px,-C1(T),x)
    xl = px*C2(T)
    x = _exp(xh,xl)
    return ldexp(x,n)
end
end
