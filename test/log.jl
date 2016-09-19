@testset for llog in (Libm.Musl.log,)
    for T in (Float32, Float64)
        @test isnan(llog(T(NaN)))
        @test llog(T(Inf)) == Inf
        @test llog(T(0)) == -Inf
        s = logspace(T(-10),T(10),1000)
        @test Base.log.(s) ≈ llog.(s)
        s = linspace(T(0),realmin(T),1000)
        @test Base.log.(s) ≈ llog.(s)
        @test_throws DomainError llog(T(-1))

        @test Base.log(prevfloat(T(Inf))) == llog(prevfloat(T(Inf)))
        @test Base.log(nextfloat(T(0))) == llog(nextfloat(T(0)))
    end
end

@testset for llog1p in (Libm.Musl.log1p,)
    for T in (Float32, Float64)
        @test isnan(llog1p(T(NaN)))
        @test llog1p(T(Inf)) == Inf
        @test llog1p(T(0)) == 0
        @test llog1p(T(-1)) == -Inf
        s = logspace(T(-10),T(10),1000)
        @test Base.log1p.(s) ≈ llog1p.(s)
        s = -logspace(T(-10),T(0),1000)
        @test Base.log1p.(s) ≈ llog1p.(s)
        @test_throws DomainError llog1p(T(-2))
        @test Base.log1p(prevfloat(T(Inf))) == llog1p(prevfloat(T(Inf)))
        @test Base.log1p(nextfloat(T(-1))) == llog1p(nextfloat(T(-1)))
    end
end
