Test = Amal
Ref = Base

@testset "denormal/nonnumber $T" for T in (Float64, Float32,)
    @test Test.exp(T(-Inf)) === T(0)
    @test Test.exp(T(Inf)) === T(Inf)
    @test Test.exp(T(NaN)) === T(NaN)
    @test Test.exp(T(0)) === T(1) # for finite argument only exp(0) is exact

    for x in T[10000, -10000]
        @test cmpdenorm(T, Test.exp(x), Ref.exp(BigFloat(x)))
    end
end

@testset "denormal/nonnumber $T" for T in (Float64, Float32,)
    @test Test.exp2(T(-Inf)) === T(0)
    @test Test.exp2(T(Inf)) === T(Inf)
    @test Test.exp2(T(NaN)) === T(NaN)
    @test Test.exp2(T(0)) === T(1) # for finite argument only exp(0) is exact

    for x in T[10000, -10000]
        @test cmpdenorm(T, Test.exp2(x), Ref.exp(BigFloat(x)))
    end
end
