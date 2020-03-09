module Libm

const FloatTypes=Union{Float32,Float64}

include("utils.jl")
include("erf.jl")
include("log/tang.jl")

export log_tang, erf, erfc


end
