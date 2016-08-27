module Libm

typealias FloatTypes Union{Float32,Float64}

include("utils.jl")
include("erf.jl")
include("log/tang.jl")

end
