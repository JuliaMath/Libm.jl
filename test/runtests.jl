using Libm
using Base.Test

isnzero{T<:AbstractFloat}(x::T) = signbit(x)
ispzero{T<:AbstractFloat}(x::T) = !signbit(x)
function cmpdenorm{Tx<:AbstractFloat, Ty<:AbstractFloat}(x::Tx, y::Ty)
    sizeof(Tx) < sizeof(Ty) ? y = Tx(y) : x = Ty(x) # cast larger type to smaller type
    (isnan(x) && isnan(y)) && return true
    (isnan(x) || isnan(y)) && return false
    (isinf(x) != isinf(y)) && return false
    (x == Tx(Inf) && y == Ty(Inf)) && return true
    (x == Tx(-Inf) && y == Ty(-Inf)) && return true
    if y == 0
        (ispzero(x) && ispzero(y)) && return true
        (isnzero(x) && isnzero(y)) && return true
        return false
    end
    (!isnan(x) && !isnan(y) && !isinf(x) && !isinf(y)) && return sign(x) == sign(y)
    return false
end

# generalize this and the following method
function ulp(x::Float64)
  x = abs(x)
  if x == 0
    return nextfloat(0.0)
  else
    exp = exponent(x)
  end
  return max(ldexp(1.0, exp-53), nextfloat(0.0))
end
function ulp(x::Float32)
  x = abs(x)
  if x == 0
    return nextfloat(0.0)
  else
    exp = exponent(x)
  end
  return max(ldexp(1.0, exp-24), nextfloat(0.0))
end

# the following compares x with y
FloatTypes = Union{Float32,Float64}
countulp{FT<:FloatTypes}(x::FT, y::BigFloat) = countulp(FT, promote(x, y)...)
function countulp(T, x::BigFloat, y::BigFloat) # wrt to y
    fx = T(x)
    fy = T(y)
    (isnan(fx) && isnan(fy)) && return 0
    (isnan(fx) || isnan(fy)) && return 10000
    if isinf(fx)
        if sign(fx) == sign(fy) && abs(fy) > 1e+300
            return 0 #Relaxed infinity handling
        end
        return 10001
    end
    (fx ===  Inf && fy ===  Inf) && return 0
    (fx === -Inf && fy === -Inf) && return 0
    if fy == 0
        if fx == 0
            return 0
        end
    return 10002
    end
    if !isnan(fx) && !isnan(fy) && !isinf(fx) && !isinf(fy)
        return abs((x - y) / ulp(T(y)))
    end
    return 10003
end

# overide domain checking that base adheres to
using Base.MPFR.ROUNDING_MODE
for f in (:sin,:cos,:tan,:asin,:acos,:tanh,:log,:asinh,:acosh,:atanh,:log10,:log1p)
    @eval begin
        import Base.$f
        function ($f)(x::BigFloat)
            z = BigFloat()
            ccall($(string(:mpfr_,f), :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32), &z, &x, ROUNDING_MODE[])
            return z
        end
    end
end


@testset "sleef" begin
for T in (Float64,) # for now only guarantee for Float64

    @testset "denormal/nonnumber" begin
        @testset "denormal/nonnumber atan2" begin
            @test xatan2(T(0.0), T(-0.0)) === T(pi)
            @test xatan2(T(-0.0), T(-0.0)) === -T(pi)
            @test ispzero(xatan2(T(0.0), T(0.0)))
            @test isnzero(xatan2(T(-0.0), T(0.0)))
            @test xatan2(T(Inf), T(-Inf)) === T(3*pi/4)
            @test xatan2(T(-Inf), T(-Inf)) === T(-3*pi/4)
            @test xatan2(T(Inf), T(Inf)) === T(pi/4)
            @test xatan2(T(-Inf), T(Inf)) === T(-pi/4)
            let
                y = T(0.0)
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
                for x in xa
                    @test xatan2(y,x) === T(pi)
                end
            end
            let
                y = T(-0.0)
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
                for x in xa
                    @test xatan2(y,x) === T(-pi)
                end
            end
            let
                ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
                xa = T[T(0.0), T(-0.0)]
                for x in xa, y in ya
                    @test xatan2(y,x) === T(-pi/2)
                end
            end
            let
                ya = T[100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5]
                xa = T[T(0.0), T(-0.0)]
                for x in xa, y in ya
                    @test xatan2(y,x) === T(pi/2)
                end
            end
            let
                y = T(Inf)
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                for x in xa
                    @test xatan2(y,x) === T(pi/2)
                end
            end
            let
                y = T(-Inf)
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                for x in xa
                    @test xatan2(y,x) === T(-pi/2)
                end
            end
            let
                ya = T[0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                x = T(Inf)
                for y in ya
                    @test ispzero(xatan2(y,x))
                end
            end
            let
                ya = T[-0.5, -1.5, -2.0, -2.5, -3.0, -100000, -100000.5]
                x = T(Inf)
                for y in ya
                    @test isnzero(xatan2(y,x))
                end
            end
            let
                ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
                x = T(NaN)
                for y in ya
                    @test isnan(xatan2(y,x))
                end
            end
            let
                y = T(NaN)
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5, NaN]
                for x in xa
                    @test isnan(xatan2(y,x))
                end
            end
        end # denormal/nonumber atan2

        @testset "denormal/nonnumber pow" begin
            @test xpow(one(T),T(NaN)) === one(T)
            @test xpow(T(NaN),zero(T)) === one(T)
            @test xpow(T(-1),T(Inf)) === one(T)
            @test xpow(T(-1),T(-Inf)) === one(T)
            let
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
                ya = T[-100000.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 100000.5]
                for x in xa, y in ya
                    @test isnan(xpow(x,y))
                end
            end
            let
                x = NaN
                ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                for y in ya
                    @test isnan(xpow(x,y))
                end
            end
            let
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5, -0.0, +0.0, 0.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                y = NaN
                for x in xa
                    @test isnan(xpow(x,y))
                end
            end
            let
                x = T(0.0)
                ya = T[1, 3, 5, 7, 100001]
                for y in ya
                    @test ispzero(xpow(x,y))
                end
            end 
            let
                x = T(-0.0)
                ya = T[1, 3, 5, 7, 100001]
                for y in ya
                    @test isnzero(xpow(x,y))
                end
            end
            let
                xa = T[T(0.0), T(-0.0)]
                ya = T[0.5, 1.5, 2.0, 2.5, 4.0, 100000, 100000.5]
                for x in xa, y in ya
                    @test ispzero(xpow(x,y))
                end
            end 
            let
                xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
                y = T(-Inf)
                for x in xa
                    @test xpow(x,y) === T(Inf)
                end
            end 
            let
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                y = T(-Inf)
                for x in xa
                    @test ispzero(xpow(x,y))
                end
            end
            let
                xa = T[-0.999, -0.5, -0.0, +0.0, +0.5, +0.999]
                y = T(Inf)
                for x in xa
                    @test ispzero(xpow(x,y))
                end
            end 
            let
                xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                y = T(Inf)
                for x in xa
                    @test xpow(x,y) === T(Inf)
                end
            end 
            let
                x = T(-Inf)
                ya = T[-100001, -5, -3, -1]
                for y in ya
                    @test isnzero(xpow(x,y))
                end
            end 
            let
                x = T(-Inf)
                ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
                for y in ya
                    @test ispzero(xpow(x,y))
                end
            end
            let
                x = T(-Inf)
                ya = T[1, 3, 5, 7, 100001]
                for y in ya
                    @test xpow(x,y) === T(-Inf)
                end
            end
            let
                x = T(-Inf)
                ya = T[0.5, 1.5, 2, 2.5, 3.5, 4, 100000, 100000.5]
                for y in ya
                    @test xpow(x,y) === T(Inf)
                end
            end
            let
                x = T(Inf)
                ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
                for y in ya
                    @test ispzero(xpow(x,y))
                end
            end
            let
                x = T(Inf)
                ya = T[0.5, 1, 1.5, 2.0, 2.5, 3.0, 100000, 100000.5]
                for y in ya
                    @test xpow(x,y) === T(Inf)
                end
            end 
            let
                x = T(0.0)
                ya = T[-100001, -5, -3, -1]
                for y in ya
                    @test xpow(x,y) ===  T(Inf)
                end
            end
            let
                x = T(-0.0)
                ya = T[-100001, -5, -3, -1]
                for y in ya
                    @test xpow(x,y) ===  T(-Inf)
                end
            end
            let
                xa = T[0.0, -0.0]
                ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
                for x in xa, y in ya
                    @test xpow(x,y) === T(Inf)
                end
            end
            let
                xa = T[1000, -1000]
                ya = T[1000, 1000.5, 1001]
                for x in xa, y in ya
                    @test cmpdenorm(xpow(x,y), BigFloat(x)^BigFloat(y))
                end
            end 
        end # denormal/nonumber pow

        @testset "denormal/nonnumber sin" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xsin(x), sin(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber sin in sincos" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    q = xsincos(x)
                    @test cmpdenorm(q.x, sin(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber cos" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xcos(x), cos(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber cos in sincos" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    q = xsincos(x)
                    @test cmpdenorm(q.y, cos(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber tan" begin
            let
                xa = T[NaN, Inf, -Inf, pi/2, -pi/2]
                for x in xa
                    @test cmpdenorm(xtan(x), tan(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber asin" begin
            let
                xa = T[NaN, Inf, -Inf, 2, -2, 1, -1]
                for x in xa
                    @test cmpdenorm(xasin(x), asin(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber acos" begin
            let
                xa = T[NaN, Inf, -Inf, 2, -2, 1, -1]
                for x in xa
                    @test cmpdenorm(xacos(x), acos(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber atan" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xatan(x), atan(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber log" begin
            let
                xa = T[NaN, Inf, -Inf, 0, -1]
                for x in xa
                    @test cmpdenorm(xlog(x), log(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber exp" begin
            let
                xa = T[NaN, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xexp(x), exp(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber sinh" begin
            let
                xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xsinh(x), sinh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber cosh" begin
            let
                xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xcosh(x), cosh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber tanh" begin
            let
                xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xtanh(x), tanh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber asinh" begin
            let
                xa = T[NaN, 0.0, -0.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xasinh(x), asinh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber acosh" begin
            let
                xa = T[NaN, 0.0, -0.0, 1.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xacosh(x), acosh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber atanh" begin
            let
                xa = T[NaN, 0.0, -0.0, 1.0, -1.0, Inf, -Inf, 10000, -10000]
                for x in xa
                    @test cmpdenorm(xatanh(x), atanh(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber cbrt" begin
            let
                xa = T[NaN, Inf, -Inf, 0.0, -0.0]
                for x in xa
                    @test cmpdenorm(xcbrt(x), cbrt(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber exp2" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xexp2(x), exp2(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber exp10" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xexp10(x), exp10(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber expm1" begin
            let
                xa = T[NaN, Inf, -Inf]
                for x in xa
                    @test cmpdenorm(xexpm1(x), expm1(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber log10" begin
            let
                xa = T[NaN, Inf, -Inf, 0.0, -1.0]
                for x in xa
                    @test cmpdenorm(xlog10(x), log10(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber log1p" begin
            let
                xa = T[NaN, Inf, -Inf, 0.0, -1.0, -2.0]
                for x in xa
                    @test cmpdenorm(xlog1p(x), log1p(BigFloat(x)))
                end
            end
        end

        @testset "denormal/nonnumber ldexp" begin
            let
                for i = -10000:10000
                    a = xldexp(T(1.0),Int32(i))
                    b = ldexp(BigFloat(1.0), i)
                    @test (isfinite(b) && a == b || cmpdenorm(a,b))
                end
            end
        end

    end #denormal/nonnumber

    @testset "accuracy (max error in ulp)" begin
        @testset "accuracy sin" begin
            rmax = 0
            xa = vcat((-10:0.0002:10, -10000000:200.1:10000000)...)
            for x in xa
                # @show x
                q = xsin(x)
                c = sin(BigFloat(x))
                u = countulp(q, c)
                rmax = max(rmax, u)
                if rmax > 1000
                    @printf("q = %.20g\nc = %.20g\nd = %.20g\nulp = %g\n", q, T(c), x, T(ulp(c)))
                end
            end
            for i = 1:10000
                s = reinterpret(Float64, reinterpret(Int64, pi/4 * i) - 20)
                e = reinterpret(Float64, reinterpret(Int64, pi/4 * i) - 20)
                d = reinterpret(Float64, reinterpret(Int64, pi/4 * i) + 1)
                for x in s:d:e
                    q = xsin(x)
                    c = sin(BigFloat(x))
                    u = countulp(q, c)
                    rmax = max(rmax, u)
                end
            end
            @printf("sin : %f ... \n", rmax)
            @test rmax < 5
        end

    end #accuracy 

end #TYPES
end #SLEEF
