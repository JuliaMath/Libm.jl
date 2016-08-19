using Libm
using BenchmarkTools

x = 2*rand(10000,1)-1

t1 = @benchmark Libm.erf.($x)
t2 = @benchmark Base.erf.($x)
t3 = @benchmark Libm.erfc.($x)
t4 = @benchmark Base.erfc.($x)

println("***** Bench erf *****")
println()
println(ratio(mean(t1),mean(t2)))
println()
println("***** Details *****")
println(t1,"\n")
println(t2,"\n")


println("***** Bench erfc *****")
println()
println(ratio(mean(t3),mean(t4)))
println()
println("***** Details *****")
println(t3,"\n")
println(t4,"\n")

