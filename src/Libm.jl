module Libm

using Base.Math.@horner

export erf, erfc, exp

include("utils.jl")
include("erf.jl")
include("exp.jl")

end
