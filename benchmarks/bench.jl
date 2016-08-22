using Libm
using BenchmarkTools

# exp benchmarks
x1 = linspace(708.4, 709.7, 10000)
x2 = -linspace(708.4, 709.7, 10000)
x3 = linspace(0.5*log(2),1.5*log(2), 10000)
x4 = -linspace(0.5*log(2),1.5*log(2), 10000)
x5 = linspace(1.5*log(2),10, 10000)
x6 = -linspace(1.5*log(2),10, 10000)
x7 = linspace(2.0^-28,2.0^-27, 10000)
x8 = -linspace(2.0^-28,2.0^-27, 10000)
x9 = linspace(1.0, 10.0, 1000)
x10 = -linspace(1.0, 10.0, 1000)
x11 = 2000*(rand(1000000)-0.5)
x = union(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)

t1 = @benchmark Libm.exp.($x)
t2 = @benchmark Base.exp.($x)

println("***** Benchmark exp *****")
println()
println(ratio(mean(t1),mean(t2)))
println()
println("Details\n")
println(t1,"\n")
println(t2,"\n")
