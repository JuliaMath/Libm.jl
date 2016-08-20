using Base.Test
using Libm

libmfun(fun) = eval(:(Libm.$fun))
basefun(fun) = eval(:(Base.$fun))


query_points = [0.5,1.0, 1.3, 1.33, 2.0, 5.0, 100.0 ,256.0, 1000.0, 2000.0]


# uniry functions
for fun in [
		:sin,:cos,
		:exp,:exp2, :log, :log10, :log2,
		:sqrt, :fabs, :ceil, :floor, :trunc, :round, :rint, :nearbyint
		]

	@testset "$fun" begin
		ourfun = libmfun(fun)
		stdfun = basefun(fun)

		for val in query_points
			@test ourfun(val) == stdfun(val)
			@test ourfun(Float32(val)) == stdfun(Float32(val))

			if fun ∉ [:log,:log10,:log2]
				@test ourfun(-val) == stdfun(-val)
				@test ourfun(Float32(-val)) == stdfun(Float32(-val))
			end
		end
	end

end

query_points2 = [-query_points; query_points ]
# binary functions
for fun in [:pow, :maxnum, :minnum]

	@testset "$fun" begin
		ourfun = libmfun(fun)
		stdfun = basefun(fun)

		for val1 in query_points2
			for val2 in query_points2
			@test ourfun(val1, val2) == stdfun(val1, val2)
			@test ourfun(Float32(val1),Float32(val2)) == stdfun(Float32(val1),Float32(val2))
		end
	end
end



