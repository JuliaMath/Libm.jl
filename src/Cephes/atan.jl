let

const TPI8  = 0.414213562373095048801688724209698078569672 # tan pi/8
const T3PI8 = 2.414213562373095048801688724209698078569672 # tan 3pi/8

# Float80 coeffs (set MBIT = 0)
# const p5 = -8.6863818178092187535440e-1
# const p4 = -1.4683508633175792446076e1
# const p3 = -6.3976888655834347413154e1
# const p2 = -9.9988763777265819915721e1
# const p1 = -5.0894116899623603312185e1

# const q6 = 1.00000000000000000000e0
# const q5 = 2.2981886733594175366172e1
# const q4 = 1.4399096122250781605352e2
# const q3 = 3.6144079386152023162701e2
# const q2 = 3.9157570175111990631099e2
# const q1 = 1.5268235069887081006606e2

# Float64 coeffs
const p5 = -8.750608600031904122785e-1
const p4 = -1.615753718733365076637e1
const p3 = -7.500855792314704667340e1
const p2 = -1.228866684490136173410e2
const p1 = -6.485021904942025371773e1

const q6 = 1.000000000000000000000e0
const q5 = 2.485846490142306297962e1
const q4 = 1.650270098316988542046e2
const q3 = 4.328810604912902668951e2
const q2 = 4.853903996359136964868e2
const q1 = 1.945506571482613964425e2

const MBIT = 6.123233995736765886130e-17 # more bits

global function atan{T}(x::T)
    if x == 0
        return x
    end
    if isinf(x)
        return copysign(T(PI2), x)
    end
    sign = false
    if signbit(x)
        sign = true
        x = -x
    end

    flag = UInt8(0)
    if x > T(T3PI8)
        flag = UInt8(1)
        y = T(PI2)
        x = -T(1.0)/x
    elseif x <= T(0.66)
        y = T(0)
    else
        flag = UInt8(2)
        y = T(PI4)
        x = (x - T(1.0)) / (x + T(1.0))
    end

    z = x*x
    # P = (@horner_oftype z p1 p2 p3 p4 p5)
    # Q = (@horner_oftype z q1 q2 q3 q4 q5 q6)

    # mod: much more accurate to use Estrin's
    zz = z*z
    P = (T(p1) + z*T(p2)) + (T(p3) + z*T(p4))*zz + T(p5)*zz*zz
    Q = (T(q1) + z*T(q2)) + (T(q3) + z*T(q4))*zz + (T(q5) + z*T(q6))*zz*zz
   
    z = x*z*P/Q + x

    if flag == UInt8(2)
        z = z + T(0.5)*T(MBIT)
    elseif flag == UInt8(1)
        z = z + T(MBIT)
    end

    y = y + z
    return sign ? -y : y
end
end


# accuracy 2.2 ulp for Float32 and Float16
# const c4 =  8.05374449538e-2
# const c3 = -1.38776856032e-1
# const c2 =  1.99777106478e-1
# const c1 = -3.33329491539e-1
# function atan{T<:SFloat}(x::T)
#     if x == 0
#         return x
#     end
#     if isinf(x)
#         return copysign(T(PI2), x)
#     end
#     sign = false
#     if signbit(x)
#         sign = true
#         x = -x
#     end

#     if x > T(T3PI8)
#         y = T(PI2)
#         x = -1/x
#     elseif x > T(TPI8)
#         y = T(PI4)
#         x = (x - T(1)) / (x + T(1))
#     else
#         y = T(0)
#     end

#     z = x*x
#     y = y + ((T(c1) + z*T(c2)) + (c3 + z*T(c4))*z*z )*z*x  + x # mod: Estrin's
#     return sign ? -y : y
# end