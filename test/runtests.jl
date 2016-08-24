using Libm
using Base.Test


nzero{T<:AbstractFloat}(::Type{T}) = -one(T)*zero(T)
pzero{T<:AbstractFloat}(::Type{T}) =  one(T)*zero(T)
isnzero{T<:AbstractFloat}(x::T) = signbit(x)
ispzero{T<:AbstractFloat}(x::T) = !signbit(x)

const TestTypes = (Float64,)

@testset "Libm" begin
for T in TestTypes

@testset "denormal/nonnumber atan2" begin
    @test xatan2(pzero(T), nzero(T)) === T(pi)
    @test xatan2(nzero(T), nzero(T)) === -T(pi)
    @test ispzero(xatan2(pzero(T), pzero(T)))
    @test isnzero(xatan2(nzero(T), pzero(T)))
    @test xatan2(T(Inf), T(-Inf)) ===  T(3*pi/4)
    @test xatan2(T(-Inf), T(-Inf)) === T(-3*pi/4)
    @test xatan2(T(Inf), T(Inf)) ===  T(pi/4)
    @test xatan2(T(-Inf), T(Inf)) ===  T(-pi/4)
    let
        y = pzero(T)
        xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
        for x in xa
            @test xatan2(y,x) === T(pi)
        end
    end
    let
        y = nzero(T)
        xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
        for x in xa
            @test xatan2(y,x) === T(-pi)
        end
    end
    let
        ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
        xa = T[pzero(T), nzero(T)]
        for x in xa, y in ya
            @test xatan2(y,x) === T(-pi/2)
        end
    end
    let
        ya = T[100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5]
        xa = T[pzero(T), nzero(T)]
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
        x = pzero(T)
        ya = T[1, 3, 5, 7, 100001]
        for y in ya
            @test ispzero(xpow(x,y))
        end
    end 
    let
        x = nzero(T)
        ya = T[1, 3, 5, 7, 100001]
        for y in ya
            @test isnzero(xpow(x,y))
        end
    end
    let
        xa = T[pzero(T), nzero(T)]
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
        x = pzero(T)
        ya = T[-100001, -5, -3, -1]
        for y in ya
            @test xpow(x,y) ===  T(Inf)
        end
    end
    let
        x = nzero(T)
        ya = T[-100001, -5, -3, -1]
        for y in ya
            @test xpow(x,y) ===  T(-Inf)
        end
    end
    let
        xa = T[pzero(T), nzero(T)]
        ya = T[-100000.5, -100000, -4, -2.5, -2, -1.5, -0.5]
        for x in xa, y in ya
            @test xpow(x,y) === T(Inf)
        end
    end
    # implement
    # let
    #     xa = T[1000, -1000]
    #     ya = T[1000, 1000.5, 1001]
    #     for x in xa, y in ya
    #         @test cmpdenorm(xpow(x,y), powfr(x,y))
    #     end
    # end 
end # denormal/nonumber pow

@testset "denormal/nonnumber sin" begin

end # denormal/nonumber test sin
end 
end
 