"""
    exp(x)

Compute the base ``e`` exponential of ``x``, in other words ``e^x``.
"""
function exp end

#  Method
#    1. Argument reduction: Reduce x to an r so that |r| <= 0.5*ln2. Given x,
#       find r and integer k such that
#
#                x = k*log2(e) + r,  |r| <= 0.5*ln2.
#
#    2. Approximate exp2(r) by a polynomial on the interval [-0.5*ln2, 0.5*ln2]:
#
#           exp(x) = 1.0 + x + polynomial(x),
#                  = polynomial(x) + x + 1 (for better accuracy)
#
#    3. Scale back: exp(x) = 2^k * exp2(r)

@inline @oftype _exp{T}(x::T) = x * x *  (0.5 + x *
    (0.166666666666666768437110590639349538832902908325195 + x *
    (4.1666666666666608842550800773096852935850620269775e-2 + x *
    (8.3333333333214711785563721946346049662679433822632e-3 + x *
    (1.3888888888923911317518911090473920921795070171356e-3 + x *
    (1.98412698864463285866599484563721489394083619117737e-4 + x *
    (2.4801587237073084572010900350491624521964695304632e-5 + x *
    (2.75572439463441191473814899370875508566314238123596e-6 + x *
    (2.7557353505471365367411173269429625065640721004456e-7 + x *
    (2.5109113077584157834518379260463349922360976052005e-8 + x *
    2.08892407222294904222385909743544413208482524169085e-9)))))))))) + x + 1

@inline @oftype _exp{T<:SmallFloatTypes}(x::T) = x * x * (0.5 + x *
    (0.1666666567325592041015625 + x *
    (4.16664779186248779296875e-2 + x *
    (8.3335675299167633056640625e-3 + x *
    (1.39336870051920413970947265625e-3 + x *
    1.9742417498491704463958740234375e-4))))) + x + 1

function exp{T}(x::T)
    # reduce
    n = _trunc(round(T(LOG2E)*x)) # truncation will give us automatic Inf handling
    r = muladd(n, -LN2U(T), x)
    r = muladd(n, -LN2L(T), r)

    # compute approximation
    u = _exp(r)
    u = _ldexp(u,n)
    
    x == T(-Inf) && (u = T(0.0))
    return u
end