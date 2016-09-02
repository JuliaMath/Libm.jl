@testset "Accuracy (max error in ulp) for $T" for T in (Float64, Float32)
 println("Accuracy tests for $T")

    fun_table = Dict(xlog => log)
    xx = map(T, vcat(0.0001:0.0001:10, 0.001:0.1:10000, 2.1.^(-1000:1000)))
    tol = 3
    test_acc(T, fun_table, xx, tol)
 

end #accuracy 