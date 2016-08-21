@testset "exp tests" begin


@test isnan(_exp(NaN))
@test isinf(_exp(Inf))
@test _exp(-Inf) == 0
@test _exp(0.0) == 1

@test isinf(_exp(710.0))
@test _exp(-746.0) == 0

x = linspace(708.4, 709.7, 25)
@test _exp.(x) == exp.(x)
x = -linspace(708.4, 709.7, 25)
@test _exp.(x) == exp.(x)
x = linspace(0.5*log(2),1.5*log(2), 25)
@test _exp.(x) == exp.(x)
x = -linspace(0.5*log(2),1.5*log(2), 25)
@test _exp.(x) == exp.(x)
x = linspace(1.5*log(2),10, 25)
@test _exp.(x) == exp.(x)
x = -linspace(1.5*log(2),10, 25)
@test _exp.(x) == exp.(x)
x = linspace(2.0^-28,2.0^-27, 25)
@test _exp.(x) == exp.(x)
x = -linspace(2.0^-28,2.0^-27, 25)
@test _exp.(x) == exp.(x)
x = linspace(1.0, 10.0, 25)
@test _exp.(x) == exp.(x)
x = -linspace(1.0, 10.0, 25)
@test _exp.(x) == exp.(x)

end
