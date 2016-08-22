@testset "scalbn" begin
    @test Libm.scalbn(1.0,Int32(0)) == 0x1p0
    @test Libm.scalbn(2.0,Int32(0)) == 0x1p1 
    @test Libm.scalbn(1.0,Int32(8)) == 0x1p8
    @test Libm.scalbn(2.0,Int32(8)) == 0x1p9
    @test Libm.scalbn(-1.0,Int32(0)) == -0x1p0
    @test Libm.scalbn(-2.0,Int32(0)) == -0x1p1
    @test Libm.scalbn(-1.0,Int32(8)) == -0x1p8
    @test Libm.scalbn(-2.0,Int32(8)) == -0x1p9

    @test Libm.scalbn(1.0,Int32(1023)) == 0x1p1023 
    @test Libm.scalbn(1.0,Int32(-1022)) == 0x1p-1022

    # ch
    @test Libm.scalbn(1.0,Int32(1024)) == Inf
    @test Libm.scalbn(1.0,Int32(1024+1023)) == Inf
    @test Libm.scalbn(1.0,Int32(1024+1023+1023)) == Inf
    @test Libm.scalbn(1.0,Int32(1024+1023+1023)) == Inf

    @test Libm.scalbn(1.0,Int32(-1023)) == 0x1p-1023
    @test Libm.scalbn(1.0,Int32(-1023-1022)) ==  0x0p0
    @test Libm.scalbn(1.0,Int32(-1023-1022-1022)) == 0x0p0
    @test Libm.scalbn(1.0,Int32(-1023-1022-1022-1022)) == 0x0p0
end