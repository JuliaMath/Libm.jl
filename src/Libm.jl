include("Musl.jl")
include("Sleef.jl")

using Reexport

@reexport using Musl
@reexport using Sleef  # won't clash since all Sleef functions are prefixed with an x
 
end