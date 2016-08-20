using Base.Test
using Libm

@testset "exp2" begin
	@test Libm.exp2(0.0) == 1.0
	@test Libm.exp2(-1.0) == 0.5
	@test Libm.exp2(1.0) == 2.0
	@test Libm.exp2(2.0) == 4





end
