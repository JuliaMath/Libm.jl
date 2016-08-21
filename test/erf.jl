@testset "erf tests" begin
    @test isnan(erf.(NaN))
    @test _erf.(Inf) == 1.0
    @test _erf.(-Inf) == -1.0
    s = linspace(-0.84375,0.84375,100)
    @test erf.(s) == _erf.(s)
    s = linspace(-2e-28,2e-28,100)
    @test erf.(s) == _erf.(s)
    s = linspace(-0.84375,0.84375,100)
    @test erf.(s) == _erf.(s)
    s = linspace(0.84375, 6, 100)
    @test erf.(s) ==  _erf.(s)
    s = linspace(-6,-0.84375, 100);
    @test erf.(s) ==  _erf.(s)
end

@testset "erfc tests" begin
    @test isnan(erfc.(NaN))
    @test _erfc.(Inf) == 0.0
    @test _erfc.(-Inf) == 2.0
    s = linspace(-0.84375,0.84375,100)
    @test erfc.(s) == _erfc.(s)
    s = linspace(-2e-56,2e-56,100)
    @test erfc.(s) == _erfc.(s)
    s = linspace(-0.25,0.25,100)
    @test erfc.(s) == _erfc.(s)
    s = linspace(0.84375, 28, 100)
    @test erfc.(s) ==  _erfc.(s)
    s = linspace(-28,-0.84375, 100);
    @test erfc.(s) ==  _erfc.(s)
end
