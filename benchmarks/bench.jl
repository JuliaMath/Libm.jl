using Libm
using BenchmarkTools

xerf1 = linspace(-0.84375,0.84375,100)
xerf2 = linspace(-2e-28,2e-28,100)
xerf3 = linspace(-0.84375,0.84375,100)
xerf4 = linspace(0.84375, 6, 100)
xerf5 = linspace(-6,-0.84375, 100);
xerf = union(xerf1,xerf2,xerf3,xerf4,xerf5);

xerfc1 = linspace(-0.84375,0.84375,100)
xerfc2 = linspace(-2e-56,2e-56,100)
xerfc3 = linspace(-0.25,0.25,100)
xerfc4 = linspace(0.84375, 28, 100)
xerfc5 = linspace(-28,-0.84375, 100);
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

