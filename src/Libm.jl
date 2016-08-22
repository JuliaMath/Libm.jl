module Libm

using Base.Math.@horner

export exp

include("exp.jl")
include("utils.jl")

end
