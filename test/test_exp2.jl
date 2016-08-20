using Base.Test
using Libm

const fuzztest = true 

@testset "exp2" begin
	@test Libm.exp2(0.0) == 1.0
	@test Libm.exp2(-1.0) == 0.5
	@test Libm.exp2(1.0) == 2.0
	@test Libm.exp2(2.0) == 4
	@test Libm.exp2(1.33) == 2.5140267490436567

	@test isnan(Libm.exp2(NaN))

	@test Libm.exp2(Inf) == Inf
	@test Libm.exp2(-Inf) == 0.0


	if fuzztest
		srand(1)
		for ii in 1: 1_000
			ff = rand(Float64)
			@test Libm.exp2(ff) == Base.exp2(ff)
		end
	end

end
