@testset "erf" begin
    for T in (Float32,Float64)
        @test isnan(Libm.Musl.erf(T(NaN)))
        @test Libm.Musl.erf(T(Inf)) == 1
        @test Libm.Musl.erf(T(-Inf)) == -1
        s = linspace(T(-0.84375),T(0.84375),100)
        @test_approx_eq Base.erf.(s) Libm.Musl.erf.(s)
        s = linspace(T(-2e-28),T(2e-28),100)
        @test_approx_eq Base.erf.(s) Libm.Musl.erf.(s)
        s = linspace(T(-0.84375),T(0.84375),100)
        @test_approx_eq Base.erf.(s) Libm.Musl.erf.(s)
        s = linspace(T(0.84375), T(6), 100)
        @test_approx_eq Base.erf.(s) Libm.Musl.erf.(s)
        s = linspace(T(-6),T(-0.84375), 100);
        @test_approx_eq Base.erf.(s) Libm.Musl.erf.(s)
    end
end

@testset "erfc" begin
    for T in (Float32,Float64)
        @test isnan(Libm.Musl.erfc(T(NaN)))
        @test Libm.Musl.erfc(T(Inf)) == 0
        @test Libm.Musl.erfc(T(-Inf)) == 2
        s = linspace(T(-0.84375),T(0.84375),100)
        @test_approx_eq Base.erfc.(s) Libm.Musl.erfc.(s)
        s = linspace(T(-2e-56),T(2e-56),100)
        @test_approx_eq Base.erfc.(s) Libm.Musl.erfc.(s)
        s = linspace(T(-0.25),T(0.25),100)
        @test_approx_eq Base.erfc.(s) Libm.Musl.erfc.(s)
        s = linspace(T(0.84375), T(28), 100)
        @test_approx_eq Base.erfc.(s) Libm.Musl.erfc.(s)
        s = linspace(T(-28),T(-0.84375), 100);
        @test_approx_eq Base.erfc.(s) Libm.Musl.erfc.(s)
    end
end
