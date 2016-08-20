using Base.Test
using Libm

libmfun(fun) = eval(:(Libm.$fun))
basefun(fun) = eval(:(Base.$fun))


query_points = [0.5,1.0, 1.3, 1.33, 2.0, 5.0, 100.0 ,256.0, 1000.0, 2000.0]


#TODO: test fabs, rint, nearbyint

#@testset "round" begin
#	@test_broken Libm.round(0.5) == 0.0
#end

# uniry functions
@testset "$fun" for fun in [
	:sin,:cos,
	:exp,:exp2, :log, :log10, :log2,
	:sqrt, :ceil, :floor, :trunc
	]


	ourfun = libmfun(fun)
	stdfun = basefun(fun)

	for val in query_points
		@test ourfun(val) == stdfun(val)
		@test ourfun(Float32(val)) == stdfun(Float32(val))

		if fun ∉ [:log,:log10,:log2,:sqrt]
			@test ourfun(-val) == stdfun(-val)
			@test ourfun(Float32(-val)) == stdfun(Float32(-val))
		end
	end
end

#TODO: Test fmax and fmin

@testset "pow" begin
	@test Libm.pow(1.0,2.0) == 1.0
	@test Libm.pow(10.0,0.0) == 1.0
	@test Libm.pow(2.0,2.0) == 4.0
	@test Libm.pow(2.0,0.5) == sqrt(2.)

	@test Libm.pow(1f0,2f0) == 1f0
	@test Libm.pow(10f0,0f0) == 1f0
	@test Libm.pow(2f0,2f0) == 4f0
	@test Libm.pow(2f0,0.5f0) == sqrt(2f0)


end


