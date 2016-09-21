
MRANGE(::Type{Float64}) = 10000000
MRANGE(::Type{Float32}) = 10000
FInt(::Type{Float64}) = Int64
FInt(::Type{Float32}) = Int32

@testset "Accuracy (max error in ulp) for $T" for T in (Float32, Float64)
 println("Accuracy tests for $T")

    fun_table = Dict(xlog => log)
    xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000)))
    tol = 3
    test_acc(T, fun_table, xx, tol)

    xx = T[]
    for i = 1:10000
        s = reinterpret(T, reinterpret(FInt(T), T(pi)/4 * i) - FInt(T)(20))
        e = reinterpret(T, reinterpret(FInt(T), T(pi)/4 * i) + FInt(T)(20))
        d = s
        while d <= e 
            append!(xx, d)
            d = reinterpret(T, reinterpret(FInt(T), d) + FInt(T)(1))
        end
    end
    xx = append!(xx, -10:0.0002:10)
    xx = append!(xx, -MRANGE(T):200.1:MRANGE(T))

    fun_table = Dict(xsin => sin, xcos => cos, xtan => tan)
    tol = 4
    test_acc(T, fun_table, xx, tol)


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


end #accuracy 