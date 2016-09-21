module Musl
export erf, erfc, log, log1p

include("Musl/utils.jl")
include("Musl/log.jl")
include("Musl/erf.jl")
end