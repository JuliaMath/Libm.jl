@testset "denormal/nonnumber $T" for T in (Float32, Float64)

@testset "denormal/nonnumber $xatan2" for xatan2 in (Libm.atan2_fast, Libm.atan2)

    @test xatan2(T(0.0),  T(-0.0)) === T(pi)
    @test xatan2(T(-0.0), T(-0.0)) === -T(pi)
    @test ispzero(xatan2(T(0.0), T(0.0)))
    @test isnzero(xatan2(T(-0.0), T(0.0)))
    @test xatan2( T(Inf), -T(Inf)) === T(3*pi/4)
    @test xatan2(-T(Inf), -T(Inf)) === T(-3*pi/4)
    @test xatan2( T(Inf),  T(Inf))  === T(pi/4)
    @test xatan2(-T(Inf),  T(Inf))  === T(-pi/4)
    

    y = T(0.0)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan2(y,x) === T(pi)
    end


    y = T(-0.0)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan2(y,x) === T(-pi)
    end


    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    xa = T[T(0.0), T(-0.0)]
    for x in xa, y in ya
        @test xatan2(y,x) === T(-pi/2)
    end


    ya = T[100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5]
    xa = T[T(0.0), T(-0.0)]
    for x in xa, y in ya
        @test xatan2(y,x) === T(pi/2)
    end


    y = T(Inf)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for x in xa
        @test xatan2(y,x) === T(pi/2)
    end


    y = T(-Inf)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for x in xa
        @test xatan2(y,x) === T(-pi/2)
    end


    ya = T[0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    x = T(Inf)
    for y in ya
        @test ispzero(xatan2(y,x))
    end


    ya = T[-0.5, -1.5, -2.0, -2.5, -3.0, -100000, -100000.5]
    x = T(Inf)
    for y in ya
        @test isnzero(xatan2(y,x))
    end


    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
    x = T(NaN)
    for y in ya
        @test isnan(xatan2(y,x))
    end


    y = T(NaN)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
    for x in xa
        @test isnan(xatan2(y,x))
    end

end # denormal/nonumber atan2


@testset "denormal/nonnumber xpow" begin

    @test Libm.pow(T(1), T(NaN))   === T(1)
    @test Libm.pow( T(NaN), T(0))  === T(1)
    @test Libm.pow(-T(1), T(Inf))  === T(1)
    @test Libm.pow(-T(1), T(-Inf)) === T(1)
    

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    ya = T[-100000.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 100000.5]
    for x in xa, y in ya
        @test isnan(Libm.pow(x,y))
    end


    x = T(NaN)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test isnan(Libm.pow(x,y))
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(NaN)
    for x in xa
        @test isnan(Libm.pow(x,y))
    end


    x = T(0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test ispzero(Libm.pow(x,y))
    end
 

    x = T(-0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test isnzero(Libm.pow(x,y))
    end


    xa = T[0.0, -0.0]
    ya = T[0.5, 1.5, 2.0, 2.5, 4.0, 100000, 100000.5]
    for x in xa, y in ya
        @test ispzero(Libm.pow(x,y))
    end
 

    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(-Inf)
    for x in xa
        @test Libm.pow(x,y) === T(Inf)
    end
 

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(-Inf)
    for x in xa
        @test ispzero(Libm.pow(x,y))
    end


    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(Inf)
    for x in xa
        @test ispzero(Libm.pow(x,y))
    end
 

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(Inf)
    for x in xa
        @test Libm.pow(x,y) === T(Inf)
    end
 

    x = T(-Inf)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test isnzero(Libm.pow(x,y))
    end
 

    x = T(-Inf)
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for y in ya
        @test ispzero(Libm.pow(x,y))
    end


    x = T(-Inf)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test Libm.pow(x,y) === T(-Inf)
    end


    x = T(-Inf)
    ya = T[0.5, 1.5, 2, 2.5, 3.5, 4, 100000, 100000.5]
    for y in ya
        @test Libm.pow(x,y) === T(Inf)
    end


    x = T(Inf)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for y in ya
        @test ispzero(Libm.pow(x,y))
    end


    x = T(Inf)
    ya = T[0.5, 1, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test Libm.pow(x,y) === T(Inf)
    end
 

    x = T(0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test Libm.pow(x,y) ===  T(Inf)
    end


    x = T(-0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test Libm.pow(x,y) ===  T(-Inf)
    end


    xa = T[0.0, -0.0]
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for x in xa, y in ya
        @test Libm.pow(x,y) === T(Inf)
    end


    xa = T[1000, -1000]
    ya = T[1000, 1000.5, 1001]
    for x in xa, y in ya
        @test cmpdenorm(Libm.pow(x,y), Base.:^(BigFloat(x),BigFloat(y)))
    end

end # denormal/nonumber pow


fun_table = Dict(Libm.sin_fast => Base.sin, Libm.sin => Base.sin)
@testset "denormal/nonnumber $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


fun_table = Dict(Libm.cos_fast => Base.cos, Libm.cos => Base.cos)
@testset "denormal/nonnumber $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber sin in $xsincos"for xsincos in (Libm.sincos_fast, Libm.sincos)
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        q = xsincos(x)[1]
        @test cmpdenorm(q, Base.sin(BigFloat(x)))
    end
end


@testset "denormal/nonnumber cos in $xsincos"for xsincos in (Libm.sincos_fast, Libm.sincos)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        q = xsincos(x)[2]
        @test cmpdenorm(q, Base.cos(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xtan" for xtan in (Libm.tan_fast, Libm.tan)
    xa = T[NaN, Inf, -Inf, pi/2, -pi/2]
    for x in xa
        @test cmpdenorm(xtan(x), Base.tan(BigFloat(x)))
    end
end


fun_table = Dict(Libm.asin => Base.asin, Libm.asin_fast => Base.asin, Libm.acos => Base.acos, Libm.acos_fast => Base.acos)
@testset "denormal/nonnumber $xatrig" for (xatrig, atrig) in fun_table
    xa = T[NaN, Inf, -Inf, 2, -2, 1, -1]
    for x in xa
        @test cmpdenorm(xatrig(x), atrig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xatan" for xatan in (Libm.atan, Libm.atan_fast)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xatan(x), Base.atan(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp" begin
    xa = T[NaN, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.exp(x), Base.exp(BigFloat(x)))
    end
end


@testset "denormal/nonnumber sinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.sinh(x), Base.sinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber cosh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.cosh(x),  Base.cosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber tanh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.tanh(x),  Base.tanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber asinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.asinh(x), Base.asinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber acosh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.acosh(x), Base.acosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber atanh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, -1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(Libm.atanh(x), Base.atanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xcbrt" for xcbrt = (Libm.cbrt, Libm.cbrt_fast)
    xa = T[NaN, Inf, -Inf, 0.0, -0.0]
    for x in xa
        @test cmpdenorm(Libm.cbrt(x), Base.cbrt(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp2" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(Libm.exp2(x), Base.exp2(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp10" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(Libm.exp10(x), Base.exp10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber expm1" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(Libm.expm1(x), Base.expm1(BigFloat(x)))
    end
end



@testset "denormal/nonnumber $xlog" for xlog in (Libm.log, Libm.log_fast)
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(xlog(x), Base.log(BigFloat(x)))
    end
end


@testset "denormal/nonnumber log10" begin
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(Libm.log10(x), Base.log10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber log2" begin
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(Libm.log2(x), Base.log2(BigFloat(x)))
    end
end


@testset "denormal/nonnumber log1p" begin
    xa = T[NaN, Inf, -Inf, 0.0, -0.0, -1.0, -2.0]
    for x in xa
        @test cmpdenorm(Libm.log1p(x), Base.log1p(BigFloat(x)))
    end
end


@testset "denormal/nonnumber ldexp" begin
    for i = -10000:10000
        a = Libm.ldexp(T(1.0), i)
        b = Base.ldexp(BigFloat(1.0), i)
        @test (isfinite(b) && a == b || cmpdenorm(a,b))
    end
end


@testset "denormal/nonnumber ilog2" begin
    @test Libm.ilog2(+T(Inf)) == typemax(Int)
    @test Libm.ilog2(-T(Inf)) == typemax(Int)
    @test Libm.ilog2(+T(0.0)) == typemin(Int)
    @test Libm.ilog2(-T(0.0)) == typemin(Int)
    @test Libm.ilog2( T(NaN)) == typemax(Int)
end


end #denormal/nonnumber
