
MRANGE(::Type{Float64}) = 10000000
MRANGE(::Type{Float32}) = 10000
IntF(::Type{Float64}) = Int64
IntF(::Type{Float32}) = Int32

@testset "Accuracy (max error in ulp) for $T" for T in (Float32, Float64)
    println("Accuracy tests for $T")

    xx = map(T, vcat(-1:0.00002:1))
    fun_table = Dict(xasin => asin, xacos => acos)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xasin_u1 => asin, xacos_u1 => acos)
    # tol = 1
    # test_acc(T, fun_table, xx, tol)

    xx = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
    fun_table = Dict(xatan => atan)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xatan_u1 => atan)
    # tol = 1
    # test_acc(T, fun_table, xx, tol)

    # xx1 = map(Tuple{T,T}, [zip(-10:0.05:10, -10:0.05:10)...])
    # xx2 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.51:100)...])
    # xx3 = map(Tuple{T,T}, [zip(-10:0.051:10, -10:0.052:10)...])
    # xx4 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.52:100)...])
    # xx = vcat(xx1, xx2, xx3, xx4)

    # fun_table = Dict(xatan2 => atan2)
    # tol = 2.5
    # test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xatan2_u1 => atan2)
    # tol = 1
    # test_acc(T, fun_table, xx, tol)

    xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    fun_table = Dict(xlog => log)
    tol = 3
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xlog_u1 => log)
    tol = 1
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xlog10 => log10)
    xx = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xlog1p => log1p)
    xx = map(T, vcat(0.0001:0.0001:10, 0.0001:0.1:10000, 10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300)))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    xx1 = map(Tuple{T,T}, [zip(-10:0.050:10, -10:0.050:10)...])
    xx2 = map(Tuple{T,T}, [zip(-10:0.051:10, -10:0.052:10)...])
    xx3 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.51:100)...])
    xx4 = map(Tuple{T,T}, [zip(-100:0.51:100, -100:0.52:100)...])
    xx = vcat(xx1, xx2, xx3, xx4)

    fun_table = Dict(xatan2 => atan2)
    tol = 2.5
    test_acc(T, fun_table, xx, tol)

    xx = T[]
    for i = 1:10000
        s = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) - IntF(T)(20))
        e = reinterpret(T, reinterpret(IntF(T), T(pi)/4 * i) + IntF(T)(20))
        d = s
        while d <= e 
            append!(xx, d)
            d = reinterpret(T, reinterpret(IntF(T), d) + IntF(T)(1))
        end
    end
    xx = append!(xx, -10:0.0002:10)
    xx = append!(xx, -MRANGE(T):200.1:MRANGE(T))

    fun_table = Dict(xsin => sin, xcos => cos, xtan => tan)
    tol = 4
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xsin_u1 => sin, xcos_u1 => cos)
    tol = 1
    test_acc(T, fun_table, xx, tol)

    sin_xsincos(x) = xsincos(x).hi
    cos_xsincos(x) = xsincos(x).lo
    fun_table = Dict(sin_xsincos => sin, cos_xsincos => cos)
    tol = 4
    test_acc(T, fun_table, xx, tol) 

    sin_xsincos_u1(x) = xsincos_u1(x).hi
    cos_xsincos_u1(x) = xsincos_u1(x).lo
    fun_table = Dict(sin_xsincos_u1 => sin, cos_xsincos_u1 => cos)
    tol = 1
    test_acc(T, fun_table, xx, tol) 


    fun_table = Dict(xpow => pow)
    xx1 = map(Tuple{T,T}, [(x,y) for x = -100:0.20:100, y = 0.1:0.20:100])[:]
    xx2 = map(Tuple{T,T}, [(x,y) for x = -100:0.21:100, y = 0.1:0.22:100])[:] # <1 f32
    xx3 = map(Tuple{T,T}, [(x,y) for x = 2.1, y = -1000:0.1:1000])
    xx = vcat(xx1, xx2, xx2)
    tol = 1.08 #WARN
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xcbrt => cbrt)
    xx = map(T, vcat(-10000:0.2:10000, 1.1.^(-1000:1000), 2.1.^(-1000:1000)))
    tol = 2
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xcbrt_u1 => cbrt)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp => exp)
    xx = map(T, vcat(-10:0.0002:10, -1000:0.1:1000))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp2 => exp2)
    xx = map(T, vcat(-10:0.0002:10, -120:0.023:1000, -1000:0.02:2000))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexp10 => exp10)
    xx = map(T, vcat(-10:0.0002:10, -35:0.023:1000, -300:0.01:300))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xexpm1 => expm1)
    xx = map(T, vcat(-10:0.0002:10, -1000:0.021:1000, -1000:0.023:1000, 
        10.0.^-(0:0.02:300), -10.0.^-(0:0.02:300), 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300)))
    tol = 2
    test_acc(T, fun_table, xx, tol)

end






    
    # fun_table = Dict(xsinh => sinh, xcosh => cosh, xtanh => tanh)
    # xx = map(T, vcat(-10:0.0002:10, -1000:0.02:1000))
    # tol = 1
    # test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xasin => asin, xacos => acos)
    # xx = map(T, vcat(-1:0.00002:1))
    # tol = 3
    # test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xexp => exp)
    # xx = map(T, vcat(-10:0.0002:10, -1000:0.1:1000))
    # tol = 1
    # test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xlog => log)
    # xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000)))
    # tol = 3
    # test_acc(T, fun_table, xx, tol)








    # fun_table = Dict(xatan => atan)
    # xx = vcat(-10:0.0002:10, -10000:0.2:10000) 
    # tol = 2.5
    # test_acc(T, fun_table, xx, tol)
    
    # xx = T[]
    # for i = 1:10000
    #     s = reinterpret(Float64, reinterpret(Int64, pi/4 * i) - 20)
    #     e = reinterpret(Float64, reinterpret(Int64, pi/4 * i) + 20)
    #     d = s
    #     while d <= e 
    #         append!(xx, d)
    #         d = reinterpret(Float64, reinterpret(Int64, d) + 1)
    #     end
    # end
    # xx = append!(xx, -10:0.0002:10)
    # xx = append!(xx, -10000000:200.1:10000000)


    # fun_table = Dict(xsin_u1 => sin, xcos_u1 => cos, xtan_u1 => tan)
    # tol = 1
    # test_acc(T, fun_table, xx, tol)
    
    # fun_table = Dict(xcbrt => cbrt)
    # xx = vcat(-10000:0.2:10000, 2.1.^(-1000:1000))
    # tol = 2
    # test_acc(T, fun_table, xx, tol)

    # fun_table = Dict(xcbrt_u1 => cbrt)
    # tol = 1
    # test_acc(T, fun_table, xx, tol)


    # fun_table = Dict(xatan2 => atan2)
    # xx1 = [(y,x) for y = -10:0.05:10, x = -10:0.05:10][:]
    # xx2 = [(y,x) for y = -100:0.51:100, x = -100:0.51:100][:]
    # xx = vcat(xx1, xx2)
    # tol = 2.5
    # test_acc(T, fun_table, xx, tol)


# end #accuracy 