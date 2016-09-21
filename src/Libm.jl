module Libm
using Reexport

include("Musl.jl")
include("Sleef.jl")

using .Musl # don't export since symbols clash
@reexport using .Sleef  # won't clash since all Sleef functions are prefixed with an x

end
