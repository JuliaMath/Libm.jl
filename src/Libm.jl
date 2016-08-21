module Libm

using Base.Math.@horner

export _erf, _erfc, _exp

include("utils.jl")
include("erf.jl")
include("exp.jl")

end
