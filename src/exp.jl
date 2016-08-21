let
global _exp

const ln2hi  =  6.93147180369123816490e-01  #0x3fe62e42, 0xfee00000 
const ln2lo  =  1.90821492927058770002e-10  #0x3dea39ef, 0x35793c76 
const invln2 =  1.44269504088896338700e+00  #0x3ff71547, 0x652b82fe
const P1     =  1.66666666666666019037e-01  #0x3FC55555, 0x5555553E 
const P2     = -2.77777777770155933842e-03  #0xBF66C16C, 0x16BEBD93 
const P3     =  6.61375632143793436117e-05  #0x3F11566A, 0xAF25DE2C 
const P4     = -1.65339022054652515390e-06  #0xBEBBBD41, 0xC5D26BF1 
const P5     =  4.13813679705723846039e-08  #0x3E663769, 0x72BEA4D0 

function _exp(x::Float64)
    hx = get_high_word(x)
    sign = hx>>31
    hx &= 0x7fffffff  # high word of |x|

    # special cases
    if hx >= 0x4086232b  # if |x| >= 708.39...
        if isnan(x)
            return x
        end
        if x > 709.782712893383973096 # overflow if x!=inf 
            x *= 0x1p1023
            return x
        end
        if x < -745.13321910194110842 # underflow if x!=-inf
            return 0.0
        end
    end

    # argument reduction
    if hx > 0x3fd62e42 # if |x| > 0.5 ln2
        if hx >= 0x3ff0a2b2  # if |x| >= 1.5 ln2
            k = unsafe_trunc(Int32, invln2*x + 0.5 - sign)
        else
            k = Int32(1) - sign - sign
        end
        hi = x - k*ln2hi  # k*ln2hi is exact here
        lo = k*ln2lo
        x = hi - lo
    elseif hx > 0x3e300000  # if |x| > 2**-28
        k = Int32(0)
        hi = x
        lo = 0
    else
        # inexact if x!=0
        return 1 + x
    end

    # x is now in primary range
    xx = x*x
    c = x - xx*(P1+xx*(P2+xx*(P3+xx*(P4+xx*P5)))) # try Horner
    y = 1 + (x*c/(2-c) - lo + hi)
    if k == Int32(0)
        return y
    end
    return scalbn(y, k)
end

end