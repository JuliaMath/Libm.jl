using Libm
using Base.Test


nzero{T<:AbstractFloat}(::Type{T}) = -one(T)*zero(T)
pzero{T<:AbstractFloat}(::Type{T}) =  one(T)*zero(T)
const TestTypes = (Float64,)

@testset "Libm" begin
# include("xexp.jl")

@testset "Denormal/nonnumber test atan2(y, x)" begin

for T in TestTypes

@test xatan2(pzero(T), nzero(T)) == T(pi)
@test xatan2(nzero(T), nzero(T)) == -T(pi)
@test !signbit(xatan2(pzero(T), pzero(T)))
@test signbit(xatan2(nzero(T), pzero(T)))
@test xatan2(T(Inf), T(-Inf)) == T(3*pi/4)
@test xatan2(T(-Inf), T(-Inf)) == T(-3*pi/4)
@test xatan2(T(Inf), T(Inf)) == T(pi/4)
@test xatan2(T(-Inf), T(Inf)) == T(-pi/4)

let
    y = pzero(T)
    xa = [-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan2(y,x) == T(pi)
    end
end

let
    y = nzero(T)
    xa = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    for x in xa
        @test xatan2(y,x) == T(-pi)
    end
end

let
    ya = T[-100000.5, -100000, -3, -2.5, -2, -1.5, -1.0, -0.5]
    xa = [pzero(T), nzero(T)]
    for x in xa, y in ya
        @test xatan2(y,x) == T(-pi/2)
    end
end

let
    ya = T[100000.5, 100000, 3, 2.5, 2, 1.5, 1.0, 0.5]
    xa = [pzero(T), nzero(T)]
    for x in xa, y in ya
        @test xatan2(y,x) == T(pi/2)
    end
end

end

end





end
 
