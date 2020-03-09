import SpecialFunctions

@testset "erf" begin
    @testset "$T" for T in (Float32,Float64)
        @test isnan(erf(T(NaN)))
        @test erf(T(Inf)) == 1
        @test erf(T(-Inf)) == -1
        s = range(T(-0.84375), T(0.84375), length=100)
        @test all(SpecialFunctions.erf.(s) .≈ erf.(s))
        s = range(T(-2e-28), T(2e-28), length=100)
        @test all(SpecialFunctions.erf.(s) .≈ erf.(s))
        s = range(T(-0.84375), T(0.84375), length=100)
        @test all(SpecialFunctions.erf.(s) .≈ erf.(s))
        s = range(T(0.84375), T(6), length=100)
        @test all(SpecialFunctions.erf.(s) .≈ erf.(s))
        s = range(T(-6), T(-0.84375), length=100)
        @test all(SpecialFunctions.erf.(s) .≈ erf.(s))
    end
end

@testset "erfc" begin
    @testset "$T" for T in (Float32,Float64)
        @test isnan(erfc(T(NaN)))
        @test erfc(T(Inf)) == 0
        @test erfc(T(-Inf)) == 2
        s = range(T(-0.84375), T(0.84375), length=100)
        @test all(SpecialFunctions.erfc.(s) .≈ erfc.(s))
        s = range(T(-2e-56), T(2e-56), length=100)
        @test all(SpecialFunctions.erfc.(s) .≈ erfc.(s))
        s = range(T(-0.25), T(0.25), length=100)
        @test all(SpecialFunctions.erfc.(s) .≈ erfc.(s))
        s = range(T(0.84375), T(28), length=100)
        @test all(SpecialFunctions.erfc.(s) .≈ erfc.(s))
        s = range(T(-28), T(-0.84375), length=100);
        @test all(SpecialFunctions.erfc.(s) .≈ erfc.(s))
    end
end
