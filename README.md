# Libm.jl

[![Travis Build Status](https://travis-ci.org/JuliaMath/Libm.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Libm.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/307l6b799amrpvks/branch/master?svg=true)](https://ci.appveyor.com/project/simonbyrne/libm-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Libm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Libm.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Libm.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Libm.jl?branch=master)

Pure Julia math library (aka libm); implementation of the functions provided by the C math library in Julia.

# Usage

We recommend running julia with `-O3` for maximal performance using `Libm.jl` and to also build a custom system image by running
```julia
# Pkg.add("WinRPM"); WinRPM.install("gcc")  # on Windows please first run this line
julia> include(joinpath(dirname(JULIA_HOME),"share","julia","build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then to restart `julia`; this will ensure you are taking full advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set)  if your CPU supports it.


To use  `Libm.jl`
```julia
julia> Pkg.clone("https://github.com/JuliaMath/Libm.jl.git")
```

Most of the useful math functions are accessed with the prefix `x`, for example:
```julia
julia> using Libm

julia> xsin(2.3)
0.7457052121767203

julia> xsin(2.3f0)
0.74570525f0

julia> xexp(3.0)
20.085536923187668

julia> xexp(3f0)
20.085537f0
```

The exported functions include (within 1 ulp)
```julia
xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
	xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp
 ```
 Faster variants include (within 4 ulp)

 ```julia

 xatan2_fast, xasin_fast, xacos_fast, xatan_fast, xsin_fast, xcos_fast, xsincos_fast,
 	 xtan_fast, xcbrt_fast, xlog_fast
```

You can also access `Libm.Musl.log(x)`  for a different implementation of the logarithmic function and `Libm.Musl.erf(x)` and `Libm.Musl.erfc(x)` for the error function and the complementary error function. 

# Benchmarks

To benchmark performance on your machine run (this will take some time and make sure you are not running anything else cpu intensive)
```julia
include(joinpath(Pkg.dir("Libm"), "bench", "bench.jl"))
```
Please run these benchmarks under `-O3` and feel free to post these results to https://github.com/JuliaMath/Libm.jl/issues/34

Sample benchmark results (commit #2c093f6)
```julia
julia> versioninfo()
Julia Version 0.5.0
Commit 3c9d753 (2016-09-19 18:14 UTC)
Platform Info:
  System: NT (x86_64-w64-mingw32)
  CPU: Intel(R) Core(TM) i7-4510U CPU @ 2.00GHz
  WORD_SIZE: 64
  BLAS: libopenblas (USE64BITINT DYNAMIC_ARCH NO_AFFINITY Haswell)
  LAPACK: libopenblas64_
  LIBM: libopenlibm
  LLVM: libLLVM-3.7.1 (ORCJIT, haswell)

sin Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.8104167251045074
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sin Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.1793971948276605
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cos Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.8349934205562576
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cos Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.1378060635629166
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tan Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.6853226996658092
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tan Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.107774919320178
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.7593163196925772
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.6007869432841224
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp2 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             4.512280527385739
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp2 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             6.361883488233164
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp10 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.31763771534986596
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp10 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.42035869376603674
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


expm1 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.2549431703992986
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


expm1 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.5746043666547447
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.9674084478517264
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.24023173518787
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log10 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.5597826044503948
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log10 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.030692929289373
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log1p Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.4666741411299218
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log1p Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.5541235936575581
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asin Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.7832199043704233
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asin Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.6338032650935848
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acos Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.557872748347827
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acos Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             4.883546516473442
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atan Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.328650262035161
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atan Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.8257447073898105
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sinh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.0352839742392277
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sinh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.0190504435521524
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cosh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.2253754945782478
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cosh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.4271157001290662
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tanh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.100372217687127
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tanh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.100830105222955
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asinh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.244364631024939
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asinh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.4517535054819344
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acosh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.21444831331582
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acosh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.2773855787139985
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atanh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.0938780085613176
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atanh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.4749757642515444
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cbrt Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.750342228215642
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cbrt Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             4.841540899764432
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


```

