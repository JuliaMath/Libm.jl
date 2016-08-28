# functions that return NaN on non-NaN argument for domain error
const libm = Base.libm_name
for f in (:atanh, :log1p)
    @eval begin
        import Base: $f
        @nowarn function ($f)(x::Float64) ccall(($(string(f)),libm), Float64, (Float64,), x) end
        @nowarn function ($f)(x::Float32) ccall(($(string(f,"f")),libm), Float32, (Float32,), x) end
    end
end
import Base.^
@nowarn function pow(x::Float64, y::Float64) ccall((:pow,libm),  Float64, (Float64,Float64), x, y) end
@nowarn function pow(x::Float32, y::Float32) ccall((:powf,libm), Float32, (Float32,Float32), x, y) end

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

    fun_table = Dict(sin => sin, cos => cos, tan => tan)
    tol = 4
    test_acc(T, fun_table, xx, tol)

 
    fun_table = Dict(asin => asin, acos => acos)
    xx = vcat(-1:0.00002:1) 
    tol = 2.5
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(atan => atan)
    xx = vcat(-10:0.0002:10, -10000:0.2:10000) 
    tol = 2.5
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(atan2 => atan2)
    xx1 = [(y,x) for y = -10:0.05:10, x = -10:0.05:10][:]
    xx2 = [(y,x) for y = -100:0.51:100, x = -100:0.51:100][:]
    xx = vcat(xx1, xx2)
    tol = 2.5
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(log => log)
    xx = vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000))
    tol = 3
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(exp => exp)
    xx = vcat(-10:0.0002:10, -1000:0.1:1000)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(cbrt => cbrt)
    xx = vcat(-10000:0.2:10000, 2.1.^(-1000:1000))
    tol = 2
    test_acc(T, fun_table, xx, tol)



    fun_table = Dict(exp2 => exp2)
    xx = vcat(-10:0.0002:10, -1000:0.02:2000)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(exp10 => exp10)
    xx = vcat(-10:0.0002:10, -300:0.01:300)
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(expm1 => expm1)
    xx = vcat(-10:0.0002:10, -1000:0.021:1000, 10.0.^(0:0.021:300), -10.0.^-(0:0.021:300))
    tol = 2
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(log10 => log10)
    xx = vcat(0.0001:0.0001:10, 0.0001:0.1:10000) 
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(log1p => log1p)
    xx = vcat(0.0001:0.0001:10, 0.0001:0.1:10000, 10.0.^(-0:0.02:300), -10.0.^(-0:0.02:300))
    tol = 1
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(sinh => sinh, cosh => cosh, tanh => tanh)
    xx = vcat(-10:0.0002:10, -1000:0.02:1000) 
    tol = 2.5
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(asinh => asinh, atanh => atanh)
    xx = vcat(-10:0.0002:10, -1000:0.02:1000) 
    tol = 2
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(acosh => acosh)
    xx = vcat(1:0.0002:10, 1:0.02:1000) 
    tol = 2
    test_acc(T, fun_table, xx, tol)


    fun_table = Dict(pow => pow)
    xx1 = [(x,y) for x = -100:0.2:100, y = 0.1:0.2:100][:]
    xx2 = [(x,y) for x = 2.1, y = -1000:0.1:1000][:]
    xx = vcat(xx1, xx2)
    tol = 1
    test_acc(T, fun_table, xx, tol)

end #accuracy 