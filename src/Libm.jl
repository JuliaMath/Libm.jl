module Libm

# package code goes here
include("erf.jl")



function GET_HIGH_WORD(d::Float64)
    bd = bits(d)
    hi = bd[1:32]
    hi = parse(UInt32,hi,2)
    return hi
end

function SET_LOW_WORD(d::Float64,lo::UInt)
    db = bits(d)
    lo = convert(UInt32,lo)
    lb = bits(lo)
    d = db[1:32]*lb
    d = parse(Int,d,2)
    d = reinterpret(Float64,d)
    return d
end


end # module
