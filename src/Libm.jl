module Libm

using Base: sign_mask, exponent_mask, significand_mask, exponent_one, exponent_bias, significand_bits

using Base.Math: @horner

export erf, erfc

include("utils.jl")
include("erf.jl")


include("log/tang.jl")

end
