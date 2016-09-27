# Libm.jl

[![Travis Build Status](https://travis-ci.org/JuliaMath/Libm.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Libm.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/307l6b799amrpvks/branch/master?svg=true)](https://ci.appveyor.com/project/simonbyrne/libm-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Libm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Libm.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Libm.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Libm.jl?branch=master)

A pure Julia math library (aka libm); implementation of the functions provided by the C math library in Julia.

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
xsin, xcos, xtan, xasin, xacos, xatan, xatan2, xsincos, xsinh, xcosh, xtanh,
    xasinh, xacosh, xatanh, xlog, xlog10, xlog1p, xilogb, xexp, xexp2, xexp10, xexpm1, xldexp, xcbrt, xpow
 ```
 Faster variants include (within 4 ulp)

 ```julia
xsin_fast, xcos_fast, xtan_fast, xsincos_fast, xasin_fast, xacos_fast, xatan_fast,
    xatan2_fast, xlog_fast, xcbrt_fast
```

You can also access `Libm.Musl.log(x)`  for a different implementation of the logarithmic function and `Libm.Musl.erf(x)` and `Libm.Musl.erfc(x)` for the error function and the complementary error function. 

# Benchmarks

You can benchmark the performance of the `Libm.jl` math library on your machine by running
```julia
include(joinpath(Pkg.dir("Libm"), "bench", "bench.jl"))
```
Note this will take some time to run for the first time (subsequent runs will be much faster). Please ensure you run the benchmarks under `-O3` and that you are not concurrently running other cpu intensive tasks.
**Feel free to post your benchmark results at https://github.com/JuliaMath/Libm.jl/issues/34**

Sample benchmark results (commit 8f94cd59f889f6a2b3e3261b968b08befb02c27e)
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
  time:             0.724478749152101
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sin Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.0660158150716432
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cos Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.8298133873989108
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cos Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.965634706162522
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tan Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.5934151370641908
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tan Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.5630016222797585
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asin Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.974368513869543
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asin Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.7409145092624554
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acos Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.0264286752944596
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acos Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.8988122562997662
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atan Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.894613021513685
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atan Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.9915863839607135
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.730954081954134
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.9217641503739673
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp2 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.5466096642613647
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp2 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             6.151655164381461
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp10 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.2731525331576639
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


exp10 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.42726460397916943
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


expm1 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.2397109747611195
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


expm1 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.4715335767109927
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.107198508628489
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.557703574314746
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log10 Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.6103944827840706
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log10 Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.03563702141266
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log1p Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.637500045965252
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


log1p Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.8174605250823508
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sinh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.9475434340518644
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


sinh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.9066813703625617
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cosh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.2877508008507996
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cosh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.3711669042945702
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tanh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.061682173979117
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


tanh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.157153334024712
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asinh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.1062075995080625
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


asinh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             1.3047918978255537
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acosh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.9450953383642957
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


acosh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             0.985581942322126
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atanh Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.2603209470598804
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


atanh Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             2.526866168660497
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cbrt Float64 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             3.125712719213493
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%


cbrt Float32 benchmark
median ratio Sleef/Base
BenchmarkTools.TrialRatio:
  time:             4.947953664071938
  gctime:           1.0
  memory:           1.0
  allocs:           1.0
  time tolerance:   5.00%
  memory tolerance: 1.00%

```

