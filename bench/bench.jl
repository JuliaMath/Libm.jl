using Libm
using JLD
using BenchmarkTools
RETUNE = false
VERBOSE = false
DETAILS = false

const bench = ("Base","Sleef","Sleef_u1")
const suite = BenchmarkGroup()
for n in bench
    suite[n] = BenchmarkGroup([n])
end

xx_sml  = logspace(-16,0,2_000_000)
xx_exp  = vcat(-20:0.00005:20, -1000:0.1:1000)
xx_log  = vcat(0.00001:0.00001:20, 0.0001:0.1:20_000, 2.0.^(-1000:1000))
xx_trig = vcat(-10:0.0000125:10, -40:0.00005:40, -40_000_000:25.0125:40_000_000)
xx_hyp  = linspace(-15, 15, 50_000_000)
xx_cbrt = vcat(-10000:0.2:10000, 2.1.^(-1000:1000))

const micros = Dict(
    "exp"   => xx_exp,
    "exp2"  => xx_exp,
    "exp10" => xx_exp,
    "expm1" => xx_sml,
    "log"   => xx_log,
    "log10" => xx_log,
    "log1p" => xx_sml,
    "sin"   => xx_trig,
    "cos"   => xx_trig,
    "tan"   => xx_trig,
    "cbrt"  => xx_cbrt,
    "sinh"  => xx_hyp,
    "cosh"  => xx_hyp,
    "tanh"  => xx_hyp
    )

const micros_u1 = Dict(
    "log"   => xx_log,
    "sin"   => xx_trig,
    "cos"   => xx_trig,
    "tan"   => xx_trig,
    "cbrt"  => xx_cbrt
    )

for n in ("Base","Sleef")
    for (f,v) in micros
        suite[n][f] = BenchmarkGroup([f])
        n == "Sleef" ? fun = Symbol("x",f) : fun = Symbol(f)
        suite[n][f] = @benchmarkable $fun.($v)
    end
end
for (f,v) in micros_u1
    suite["Sleef_u1"][f] = BenchmarkGroup([f])
    fun = Symbol("x",f,"_u1")
    suite["Sleef_u1"][f] = @benchmarkable $fun.($v)
end

paramf = joinpath(dirname(@__FILE__), "params.jld")
if !isfile(paramf) || RETUNE
    tune!(suite; verbose=VERBOSE)
    save(paramf, "suite", params(suite))
    println("Saving tuned parameters.")
else
    println("Loading pretuned parameters.")
    loadparams!(suite, load(paramf, "suite"), :evals, :samples)
end

println("Warming up...")
warmup(suite,VERBOSE)
println("Running micro benchmarks...")
results = run(suite; verbose=VERBOSE)

for f in sort(collect(keys(micros)))
    println()
    print_with_color(:magenta, string(f, " benchmark\n"))
    print_with_color(:blue, "median ratio Sleef/Base\n")
    println(ratio(median(results["Sleef"][f]), median(results["Base"][f])))
    println()
    if DETAILS
        print_with_color(:blue, "details Sleef/Base\n")
        println(results["Sleef"][f])
        println(results["Base"][f])
        println()
    end
end
for f in sort(collect(keys(micros_u1)))
    println()
    print_with_color(:magenta, string(f, " benchmark\n"))
    print_with_color(:blue, "median ratio Sleef_u1/Base\n")
    println(ratio(median(results["Sleef_u1"][f]), median(results["Base"][f])))
    println()
    if DETAILS
        print_with_color(:blue, "details Sleef_u1/Base\n")
        println(results["Sleef_u1"][f])
        println(results["Base"][f])
        println()
    end
end