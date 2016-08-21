using Libm
using BenchmarkTools

# erf and benchmarks
xerf1 = linspace(-0.84375,0.84375,10000)
xerf2 = linspace(-2e-28,2e-28,10000)
xerf3 = linspace(-0.84375,0.84375,10000)
xerf4 = linspace(0.84375, 6,10000)
xerf5 = linspace(-6,-0.84375,10000);
xerf = union(xerf1,xerf2,xerf3,xerf4,xerf5);

xerfc1 = linspace(-0.84375,0.84375,10000)
xerfc2 = linspace(-2e-56,2e-56,10000)
xerfc3 = linspace(-0.25,0.25,10000)
xerfc4 = linspace(0.84375, 28,10000)
xerfc5 = linspace(-28,-0.84375,10000);
xercf = union(xerfc1,xerfc2,xerfc3,xerfc4,xerfc5)

t1 = @benchmark Libm.erf.($xerf)
t2 = @benchmark Base.erf.($xerf)
t3 = @benchmark Libm.erfc.($xercf)
t4 = @benchmark Base.erfc.($xercf)

println("***** Benchmark erf *****")
println()
println(ratio(mean(t1),mean(t2)))
println()
println("Details\n")
println(t1,"\n")
println(t2,"\n")


println("***** Benchmark erfc *****")
println()
println(ratio(mean(t3),mean(t4)))
println()
println("Details\n")
println(t3,"\n")
println(t4,"\n")

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

xexp = union(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);

t1 = @benchmark Libm._exp.($xexp)
t2 = @benchmark Base.exp.($xexp)

println("***** Benchmark exp *****")
println()
println(ratio(mean(t1),mean(t2)))
println()
println("Details\n")
println(t1,"\n")
println(t2,"\n")
