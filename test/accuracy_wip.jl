@testset "Accuracy (max error in ulp) for $T" for T in (Float64, )
 println("Accuracy tests for $T")

    # fun_table = Dict(xlog => log)
    # xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000)))
    # tol = 3
    # test_acc(T, fun_table, xx, tol)


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


    fun_table = Dict(xsin_u1 => sin, xcos_u1 => cos, xtan_u1 => tan)
    tol = 1
    test_acc(T, fun_table, xx, tol)
    
    fun_table = Dict(xcbrt => cbrt)
    xx = vcat(-10000:0.2:10000, 2.1.^(-1000:1000))
    tol = 2
    test_acc(T, fun_table, xx, tol)

    fun_table = Dict(xcbrt_u1 => cbrt)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(xatan2 => atan2)
    xx1 = [(y,x) for y = -10:0.05:10, x = -10:0.05:10][:]
    xx2 = [(y,x) for y = -100:0.51:100, x = -100:0.51:100][:]
    xx = vcat(xx1, xx2)
    tol = 2.5
    test_acc(T, fun_table, xx, tol)


end #accuracy 