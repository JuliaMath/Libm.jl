@testset "denormal/nonnumber $T" for T in (Float32, Float64)

@testset "denormal/nonnumber $xatan2" for xatan2 in (xatan2_fast, xatan2)

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

    @test xpow(T(1), T(NaN))   === T(1)
    @test xpow( T(NaN), T(0))  === T(1)
    @test xpow(-T(1), T(Inf))  === T(1)
    @test xpow(-T(1), T(-Inf)) === T(1)
    

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    ya = T[-100000.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 100000.5]
    for x in xa, y in ya
        @test isnan(xpow(x,y))
    end


    x = T(NaN)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test isnan(xpow(x,y))
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(NaN)
    for x in xa
        @test isnan(xpow(x,y))
    end


    x = T(0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test ispzero(xpow(x,y))
    end
 

    x = T(-0.0)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test isnzero(xpow(x,y))
    end


    xa = T[0.0, -0.0]
    ya = T[0.5, 1.5, 2.0, 2.5, 4.0, 100000, 100000.5]
    for x in xa, y in ya
        @test ispzero(xpow(x,y))
    end
 

    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(-Inf)
    for x in xa
        @test xpow(x,y) === T(Inf)
    end
 

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(-Inf)
    for x in xa
        @test ispzero(xpow(x,y))
    end


    xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
    y = T(Inf)
    for x in xa
        @test ispzero(xpow(x,y))
    end
 

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = T(Inf)
    for x in xa
        @test xpow(x,y) === T(Inf)
    end
 

    x = T(-Inf)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test isnzero(xpow(x,y))
    end
 

    x = T(-Inf)
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for y in ya
        @test ispzero(xpow(x,y))
    end


    x = T(-Inf)
    ya = T[1, 3, 5, 7, 100001]
    for y in ya
        @test xpow(x,y) === T(-Inf)
    end


    x = T(-Inf)
    ya = T[0.5, 1.5, 2, 2.5, 3.5, 4, 100000, 100000.5]
    for y in ya
        @test xpow(x,y) === T(Inf)
    end


    x = T(Inf)
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for y in ya
        @test ispzero(xpow(x,y))
    end


    x = T(Inf)
    ya = T[0.5, 1, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test xpow(x,y) === T(Inf)
    end
 

    x = T(0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test xpow(x,y) ===  T(Inf)
    end


    x = T(-0.0)
    ya = T[-100001, -5, -3, -1]
    for y in ya
        @test xpow(x,y) ===  T(-Inf)
    end


    xa = T[0.0, -0.0]
    ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
    for x in xa, y in ya
        @test xpow(x,y) === T(Inf)
    end


    xa = T[1000, -1000]
    ya = T[1000, 1000.5, 1001]
    for x in xa, y in ya
        @test cmpdenorm(xpow(x,y), pow(BigFloat(x),BigFloat(y)))
    end

end # denormal/nonumber pow


fun_table = Dict(xsin_fast => sin, xsin => sin)
@testset "denormal/nonnumber $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


fun_table = Dict(xcos_fast => cos, xcos => cos)
@testset "denormal/nonnumber $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber sin in $xsincos"for xsincos in (xsincos_fast, xsincos)
    xa = T[NaN, -0.0, 0.0, Inf, -Inf]
    for x in xa
        q = xsincos(x)
        @test cmpdenorm(q.hi, sin(BigFloat(x)))
    end
end


@testset "denormal/nonnumber cos in $xsincos"for xsincos in (xsincos_fast, xsincos)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        q = xsincos(x)
        @test cmpdenorm(q.lo, cos(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xtan" for xtan in (xtan_fast, xtan)
    xa = T[NaN, Inf, -Inf, pi/2, -pi/2]
    for x in xa
        @test cmpdenorm(xtan(x), tan(BigFloat(x)))
    end
end


fun_table = Dict(xasin => asin, xasin_fast => asin, xacos => acos, xacos_fast => acos)
@testset "denormal/nonnumber $xatrig" for (xatrig, atrig) in fun_table
    xa = T[NaN, Inf, -Inf, 2, -2, 1, -1]
    for x in xa
        @test cmpdenorm(xatrig(x), atrig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xatan" for xatan in (xatan, xatan_fast)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xatan(x), atan(BigFloat(x)))
    end
end



@testset "denormal/nonnumber $xlog" for xlog in (xlog, xlog_fast)
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(xlog(x), log(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xexp" begin
    xa = T[NaN, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xexp(x), exp(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xsinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xsinh(x), sinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xcosh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xcosh(x), cosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xtanh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xtanh(x), tanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xasinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xasinh(x), asinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xacosh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xacosh(x), acosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xatanh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, -1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xatanh(x), atanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xcbrt" for xcbrt = (xcbrt, xcbrt_fast)
    xa = T[NaN, Inf, -Inf, 0.0, -0.0]
    for x in xa
        @test cmpdenorm(xcbrt(x), cbrt(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xexp2" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexp2(x), exp2(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xexp10" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexp10(x), exp10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xexpm1" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexpm1(x), expm1(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xlog10" begin
    xa = T[NaN, Inf, -Inf, 0.0, -1.0]
    for x in xa
        @test cmpdenorm(xlog10(x), log10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xlog1p" begin
    xa = T[NaN, Inf, -Inf, 0.0, -0.0, -1.0, -2.0]
    for x in xa
        @test cmpdenorm(xlog1p(x), log1p(BigFloat(x)))
    end
end


@testset "denormal/nonnumber xldexp" begin
    for i = -10000:10000
        a = xldexp(T(1.0), i)
        b = ldexp(BigFloat(1.0), i)
        @test (isfinite(b) && a == b || cmpdenorm(a,b))
    end
end


@testset "denormal/nonnumber xilogb" begin
    @test xilogb(+T(Inf)) == typemax(Int)
    @test xilogb(-T(Inf)) == typemax(Int)
    @test xilogb(+T(0.0)) == typemin(Int)
    @test xilogb(-T(0.0)) == typemin(Int)
    @test xilogb( T(NaN)) == typemax(Int)
end


end #denormal/nonnumber
