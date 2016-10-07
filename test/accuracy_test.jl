# using Sleef
# using Cephes


const Test = Amal
const Ref = Base

@testset "Accuracy (max error in ulp) for $T" for T in (Float64 , Float32, Float16, )
    println("Accuracy tests for $T")


    xx = map(T, linspace(-3,3,10000000))
    fun_table = Dict(Test.exp => Ref.exp)
    tol = 5
    test_acc(T, fun_table, xx, tol)


    # # xx = map(T, vcat(-10:0.0002:10, -10000:0.2:10000, -10000:0.201:10000))
    # # xx = map(T, vcat(-10:0.0002:10,-10:0.0001:10, -2:0.00001:2))
    # xx = map(T, vcat(-1:0.00001:1, -1:0.000001:1,-1:0.000021:1))
    # fun_table = Dict(Test.atan => Ref.atan)
    # tol = 3
    # test_acc(T, fun_table, xx, tol)

end
