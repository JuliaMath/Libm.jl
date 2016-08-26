@testset "denormal/nonnumber $T" for T in (Float64,)

@testset "denormal/nonnumber $xatan2" for xatan2 in (xatan2, xatan2_u1)

    @test xatan2(T(0.0), T(-0.0)) === T(pi)
    @test xatan2(T(-0.0), T(-0.0)) === -T(pi)
    @test ispzero(xatan2(T(0.0), T(0.0)))
    @test isnzero(xatan2(T(-0.0), T(0.0)))
    @test xatan2(T(Inf), T(-Inf)) === T(3*pi/4)
    @test xatan2(T(-Inf), T(-Inf)) === T(-3*pi/4)
    @test xatan2(T(Inf), T(Inf)) === T(pi/4)
    @test xatan2(T(-Inf), T(Inf)) === T(-pi/4)
    

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


@testset "denormal/nonnumber pow" begin

    @test xpow(one(T),T(NaN)) === one(T)
    @test xpow(T(NaN),zero(T)) === one(T)
    @test xpow(T(-1),T(Inf)) === one(T)
    @test xpow(T(-1),T(-Inf)) === one(T)
    

    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    ya = T[-100000.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 100000.5]
    for x in xa, y in ya
        @test isnan(xpow(x,y))
    end


    x = NaN
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    for y in ya
        @test isnan(xpow(x,y))
    end


    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
    y = NaN
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


    xa = T[T(0.0), T(-0.0)]
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


fun_table = Dict(xsin => sin, xsin_u1 => sin, xcos => cos, xcos_u1 => cos)
@testset "denormal/nonnumber $xtrig" for (xtrig, trig) in fun_table
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xtrig(x), trig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $trig in $xsincos"for xsincos in (xsincos, xsincos_u1), trig in (sin, cos)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        q = xsincos(x)
        trig == sin ? val = q.x : val = q.y
        @test cmpdenorm(val, trig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xtan" for xtan in (xtan, xtan_u1)
    xa = T[NaN, Inf, -Inf, pi/2, -pi/2]
    for x in xa
        @test cmpdenorm(xtan(x), tan(BigFloat(x)))
    end
end


fun_table = Dict(xasin => asin, xasin_u1 => asin, xacos => acos, xacos_u1 => acos)
@testset "denormal/nonnumber $xatrig" for (xatrig, atrig) in fun_table
    xa = T[NaN, Inf, -Inf, 2, -2, 1, -1]
    for x in xa
        @test cmpdenorm(xatrig(x), atrig(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xatan" for xatan in (xatan, xatan_u1)
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xatan(x), atan(BigFloat(x)))
    end
end



@testset "denormal/nonnumber $xlog" for xlog in (xlog, xlog_u1)
    xa = T[NaN, Inf, -Inf, 0, -1]
    for x in xa
        @test cmpdenorm(xlog(x), log(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp" begin
    xa = T[NaN, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xexp(x), exp(BigFloat(x)))
    end
end


@testset "denormal/nonnumber sinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xsinh(x), sinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber cosh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xcosh(x), cosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber tanh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xtanh(x), tanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber asinh" begin
    xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xasinh(x), asinh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber acosh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xacosh(x), acosh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber atanh" begin
    xa = T[NaN, 0.0, -0.0, 1.0, -1.0, Inf, -Inf, 10000, -10000]
    for x in xa
        @test cmpdenorm(xatanh(x), atanh(BigFloat(x)))
    end
end


@testset "denormal/nonnumber $xcbrt" for xcbrt = (xcbrt, xcbrt_u1)
    xa = T[NaN, Inf, -Inf, 0.0, -0.0]
    for x in xa
        @test cmpdenorm(xcbrt(x), cbrt(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp2" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexp2(x), exp2(BigFloat(x)))
    end
end


@testset "denormal/nonnumber exp10" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexp10(x), exp10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber expm1" begin
    xa = T[NaN, Inf, -Inf]
    for x in xa
        @test cmpdenorm(xexpm1(x), expm1(BigFloat(x)))
    end
end


@testset "denormal/nonnumber log10" begin
    xa = T[NaN, Inf, -Inf, 0.0, -1.0]
    for x in xa
        @test cmpdenorm(xlog10(x), log10(BigFloat(x)))
    end
end


@testset "denormal/nonnumber log1p" begin
    xa = T[NaN, Inf, -Inf, 0.0, -1.0, -2.0]
    for x in xa
        @test cmpdenorm(xlog1p(x), log1p(BigFloat(x)))
    end
end


@testset "denormal/nonnumber ldexp" begin
    for i = -10000:10000
        a = xldexp(T(1.0),Int32(i))
        b = ldexp(BigFloat(1.0), i)
        @test (isfinite(b) && a == b || cmpdenorm(a,b))
    end
end

end #denormal/nonnumber
