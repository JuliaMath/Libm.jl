module Libm

include("Musl/Musl.jl")
include("Sleef/Sleef.jl")
include("Cephes/Cephes.jl")

using .Musl
using .Sleef
using .Cephes

end