using Libm
using BenchmarkTools

x = 2*rand(10000,1)-1

println("Bench erf")
t1 = @benchmark Libm.erf.($x)
t2 = @benchmark Base.erf.($x)
println(ratio(mean(t1),mean(t2)))

println()
println("Bench erfc")
t1 = @benchmark Libm.erfc.($x)
t2 = @benchmark Base.erfc.($x)
println(ratio(mean(t1),mean(t2)))

