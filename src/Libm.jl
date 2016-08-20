module Libm
__precompile__()

using Base.Math.@horner

include("utils.jl")
include("erf.jl")
include("exp.jl")

end
