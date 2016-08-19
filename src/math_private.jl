#Relates to OpenLibm: src/math_private, however little of that code can be transfered as it uses C unions.



function highword(d::Float64)
	(reinterpret(UInt64, d) >> 32) % UInt32

function lowword(d::Float64)
	reinterpret(UInt64, d) % UInt32



