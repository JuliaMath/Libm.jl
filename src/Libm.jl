module Libm

using Base.Math.@horner

export erf, erfc

include("utils.jl")
include("erf.jl")

end
