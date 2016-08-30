@testset "Accuracy (max error in ulp) for $T" for T in (Float64,)
println("Accuracy (max error in ulp) for $T")

    xx = T[]
    for i = 1:10000
        s = reinterpret(Float64, reinterpret(Int64, pi/4 * i) - 20)
        e = reinterpret(Float64, reinterpret(Int64, pi/4 * i) + 20)
        d = s
        while d <= e 
            append!(xx, d)
            d = reinterpret(Float64, reinterpret(Int64, d) + 1)
        end
    end
    xx = append!(xx, -10:0.0002:10)
    xx = append!(xx, -10000000:200.1:10000000)

    fun_table = Dict(xsin => sin, xcos => cos, xtan => tan)
    tol = 4
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xsin_u1 => sin, xcos_u1 => cos, xtan_u1 => tan)
    tol = 1
    test_acc(T, fun_table, xx, tol)

    sin_xsincos(x) = xsincos(x).x
    cos_xsincos(x) = xsincos(x).y
    fun_table = Dict(sin_xsincos => sin, cos_xsincos => cos)
    tol = 4
    test_acc(T, fun_table, xx, tol) 

    sin_xsincos_u1(x) = xsincos_u1(x).x
    cos_xsincos_u1(x) = xsincos_u1(x).y
    fun_table = Dict(sin_xsincos_u1 => sin, cos_xsincos_u1 => cos)
    tol = 1
    test_acc(T, fun_table, xx, tol) 


    fun_table = Dict(xasin => asin, xacos => acos)
    xx = vcat(-1:0.00002:1) 
    tol = 2.5
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xasin_u1 => asin, xacos_u1 => acos)
    xx = vcat(-1:0.00002:1) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xatan => atan)
    xx = vcat(-10:0.0002:10, -10000:0.2:10000) 
    tol = 2.5
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xatan_u1 => atan)
    xx = vcat(-10:0.0002:10, -10000:0.2:10000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xatan2 => atan2)
    xx1 = [(y,x) for y = -10:0.05:10, x = -10:0.05:10][:]
    xx2 = [(y,x) for y = -100:0.51:100, x = -100:0.51:100][:]
    xx = vcat(xx1, xx2)
    tol = 2.5
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xatan2_u1 => atan2)
    xx = vcat(xx1, xx2)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xlog => log)
    xx = vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000))
    tol = 3
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xlog_u1 => log)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp => exp)
    xx = vcat(-10:0.0002:10, -1000:0.1:1000)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xcbrt => cbrt)
    xx = vcat(-10000:0.2:10000, 2.1.^(-1000:1000))
    tol = 2
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xcbrt_u1 => cbrt)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp2 => exp2)
    xx = vcat(-10:0.0002:10, -1000:0.02:2000)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp10 => exp10)
    xx = vcat(-10:0.0002:10, -300:0.01:300)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexpm1 => expm1)
    xx = vcat(-10:0.0002:10, -1000:0.021:1000, 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300))
    tol = 2
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xlog10 => log10)
    xx = vcat(0.0001:0.0001:10, 0.0001:0.1:10000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xlog1p => log1p)
    xx = vcat(0.0001:0.0001:10, 0.0001:0.1:10000, 10.0.^(-0:0.02:300), -10.0.^(-0:0.02:300))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xsinh => sinh, xcosh => cosh, xtanh => tanh)
    xx = vcat(-10:0.0002:10, -1000:0.02:1000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xasinh => asinh, xatanh => atanh)
    xx = vcat(-10:0.0002:10, -1000:0.02:1000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xacosh => acosh)
    xx = vcat(1:0.0002:10, 1:0.02:1000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xpow => pow)
    xx1 = [(x,y) for x = -100:0.2:100, y = 0.1:0.2:100][:]
    xx2 = [(x,y) for x = 2.1, y = -1000:0.1:1000][:]
    xx = vcat(xx1, xx2)
    tol = 1
    test_acc(T, fun_table, xx, tol)

    @testset "xilogb at arbitrary values" begin
        xd = Dict{T,Int}(T(1e-30) => -100, T(2.31e-11) => -36, T(-1.0) => 0, T(1.0) => 0, 
                    T(2.31e11) => 37,  T(1e30) => 99)
        for (i,j) in xd
            @test xilogb(i)  === j
        end
    end
    
end #accuracy 