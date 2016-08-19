@testset "erf" begin
    @test isnan(Libm.erf(NaN))
    @test Libm.erf(Inf) == 1
    @test Libm.erf(-Inf) == -1
    s = linspace(-0.84375,0.84375,100)
    @test_approx_eq Base.erf.(s) Libm.erf.(s)
    s = linspace(-2e-28,2e-28,100)
    @test_approx_eq Base.erf.(s) Libm.erf.(s)
    s = linspace(-0.84375,0.84375,100)
    @test_approx_eq Base.erf.(s) Libm.erf.(s)
    s = linspace(0.84375, 6, 100)
    @test_approx_eq Base.erf.(s) Libm.erf.(s)
    s = linspace(-6,-0.84375, 100);
    @test_approx_eq Base.erf.(s) Libm.erf.(s)
end

@testset "erfc" begin
    @test isnan(Libm.erfc(NaN))
    @test Libm.erfc(Inf) == 0
    @test Libm.erfc(-Inf) == 2
    s = linspace(-0.84375,0.84375,100)
    @test_approx_eq Base.erfc.(s) Libm.erfc.(s)
    s = linspace(-2e-56,2e-56,100)
    @test_approx_eq Base.erfc.(s) Libm.erfc.(s)
    s = linspace(-0.25,0.25,100)
    @test_approx_eq Base.erfc.(s) Libm.erfc.(s)
    s = linspace(0.84375, 28, 100)
    @test_approx_eq Base.erfc.(s) Libm.erfc.(s)
    s = linspace(-28,-0.84375, 100);
    @test_approx_eq Base.erfc.(s) Libm.erfc.(s)
end
