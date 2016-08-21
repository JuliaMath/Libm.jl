module Libm

using Base.Math.@horner

export _erf, _erfc

include("utils.jl")
include("erf.jl")

end
