using Base.Math: DomainError

@testset "log" begin
  @testset "$llog" for llog in (log_tang,)
      for T in (Float32, Float64)
          @test isnan(llog(T(NaN)))
          @test llog(T(Inf)) == Inf
          @test llog(T(0)) == -Inf
          s = exp10.(range(T(-10), stop=T(10), length=1000))
          @test Base.log.(s) ≈ llog.(s)
          s = range(T(0), stop=floatmin(T), length=1000)
          @test Base.log.(s) ≈ llog.(s)
          @test_throws DomainError llog(T(-1))

          @test Base.log(prevfloat(T(Inf))) == llog(prevfloat(T(Inf)))
          @test Base.log(nextfloat(T(0))) == llog(nextfloat(T(0)))
      end
  end

  @testset "$llog1p" for llog1p in (log1p_tang,)
      for T in (Float32, Float64)
          @test isnan(llog1p(T(NaN)))
          @test llog1p(T(Inf)) == Inf
          @test llog1p(T(0)) == 0
          @test llog1p(T(-1)) == -Inf
          s = exp10.(range(T(-10), stop=T(10), length=1000))
          @test Base.log1p.(s) ≈ llog1p.(s)
          s = -exp10.(range(T(-10), stop=T(0), length=1000))
          @test Base.log1p.(s) ≈ llog1p.(s)
          @test_throws DomainError llog1p(T(-2))
          @test Base.log1p(prevfloat(T(Inf))) == llog1p(prevfloat(T(Inf)))
          @test Base.log1p(nextfloat(T(-1))) == llog1p(nextfloat(T(-1)))
      end
  end
end
