@testset for llog in (Libm.log_tang,)
    for T in (Float32, Float64)
        @test isnan(llog(T(NaN)))
        @test llog(T(Inf)) == Inf
        @test llog(T(0)) == -Inf
        s = logspace(T(-10),T(10),1000)
        @test Base.log.(s) ≈ llog.(s)
        s = linspace(T(0),realmin(T),1000)
        @test Base.log.(s) ≈ llog.(s)
        @test_throws DomainError llog(T(-1))
    end
end

@testset for llog1p in (Libm.log1p_tang,)
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
    end
end
