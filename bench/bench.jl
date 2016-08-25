using Libm
using JLD
using BenchmarkTools

const bench = ("Base","Libm")
RETUNE = false
VERBOSE = false

const suite = BenchmarkGroup()
for n in bench
    suite[n] = BenchmarkGroup([n])
end
srand(100)

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
x11 = 40000*(rand(2_000_000)-0.5)
xx_exp = union(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11)
xx_sm = logspace(-16,0,2_000_000)
xx_log = 4000*rand(2_000_000)
xx_trig =  vcat((-40:0.00005:40, -40_000_000:25.0125:40_000_000, 1:0.125:40000)...)
xx_hyp = linspace(-15, 15, 50_000_000)

const micros = Dict(
    # "exp"   => xx_exp,
    # "exp2"  => xx_exp,
    # "exp10" => xx_exp,
    # "expm1" => xx_sm,
    # "log"   => xx_log,
    # "log10" => xx_log,
    # "log1p" => xx_sm,
    "sin"   => xx_trig,
    "cos"   => xx_trig,
    "tan"   => xx_trig
    # "sinh"  => xx_hyp,
    # "cosh"  => xx_hyp,
    # "tanh"  => xx_hyp
    )


for n in keys(suite)
    for (f,v) in micros
        suite[n][f] = BenchmarkGroup([f])
        n == "Libm" ? fun = Symbol("x",f) : fun = Symbol(f)
        suite[n][f] = @benchmarkable $fun.($v)
    end
end
paramf = joinpath("./", "bench", "params.jld") #fixme
if !isfile(paramf) || RETUNE
    tune!(suite, verbose=true)
    save(paramf, "suite", params(suite))
    println("Saving tuned parameters.")
else
    println("Loading pretuned parameters.")
    loadparams!(suite, load(paramf, "suite"), :evals, :samples)
end
println("Warming up...")
warmup(suite)
println("Running micro benchmarks...")
results = run(suite, verbose=true)

for f in sort(collect(keys(micros)))
    println()
    print_with_color(:magenta, string(f, " benchmark\n"))
    print_with_color(:blue, "median ratio Libm/Base\n")
    println(ratio(median(results["Libm"][f]), median(results["Base"][f])))
    println()
    if VERBOSE
    print_with_color(:blue, "details Libm/Base\n")
    println(results["Libm"][f])
    println(results["Base"][f])
    println()
    end
end