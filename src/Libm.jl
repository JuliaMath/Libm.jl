module Libm

include("Musl/Musl.jl")
include("Sleef/Sleef.jl")
include("Cephes/Cephes.jl")
include("Amal/Amal.jl")

using .Musl
using .Sleef
using .Cephes
using .Amal

end