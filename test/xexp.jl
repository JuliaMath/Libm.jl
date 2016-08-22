@testset "xexp" begin
    @test isnan(xexp(NaN))
    @test isinf(xexp(Inf))
    @test xexp(-Inf) == 0
    @test xexp(0.0) == 1

    @test isinf(xexp(710.0))
    @test xexp(-746.0) == 0

    x = linspace(708.4, 709.7, 25)
    @test_approx_eq xexp.(x) exp.(x)

    x = -linspace(708.4, 709.7, 25)
    @test_approx_eq xexp.(x) exp.(x)

    x = linspace(0.5*log(2),1.5*log(2), 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = -linspace(0.5*log(2),1.5*log(2), 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = linspace(1.5*log(2),10, 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = -linspace(1.5*log(2),10, 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = linspace(2.0^-28,2.0^-27, 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = -linspace(2.0^-28,2.0^-27, 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = linspace(1.0, 10.0, 25)
    @test_approx_eq xexp.(x) exp.(x)
    
    x = -linspace(1.0, 10.0, 25)
    @test_approx_eq xexp.(x) exp.(x)
end