using Libm
using BenchmarkTools

suite = BenchmarkGroup(["Libm"])

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
xx_trig = linspace(-100, 100, 50_000_000)
xx_hyp = linspace(-15, 15, 50_000_000)

const micros = Dict(
    "exp"   => xx_exp,
    "exp2"  => xx_exp,
    "exp10" => xx_exp,
    "expm1" => xx_sm,
    "log"   => xx_log,
    "log10" => xx_log,
    "log1p" => xx_sm,
    "sin"   => xx_trig,
    "cos"   => xx_trig,
    "tan"   => xx_trig,
    "sinh"  => xx_hyp,
    "cosh"  => xx_hyp,
    "tanh"  => xx_hyp
    )

for (f,v) in micros
    suite[string(f)] = BenchmarkGroup([string(f)])
    libmf = Symbol("x",f)
    basef = Symbol(f)
    suite[f]["Libm"] = @benchmarkable $libmf.($v)
    suite[f]["Base"] = @benchmarkable $basef.($v)
end

tune!(suite)
results = run(suite)

for f in sort(collect(keys(micros)))
    println()
    print_with_color(:magenta, string(f, " benchmark\n"))
    print_with_color(:blue, "median ratio Libm/Base\n")
    println(ratio(median(result[f]["Libm"]), median(result[f]["Base"])))
    println()
    # print_with_color(:blue, "details Libm/Base\n")
    # println(result[f]["Libm"])
    # println(result[f]["Base"])
    # println()
end