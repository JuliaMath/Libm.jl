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

julia> Libm.sin(2.3)
0.7457052121767203

julia> Libm.sin(2.3f0)
0.74570525f0

julia> Libm.exp(3.0)
20.085536923187668

julia> Libm.exp(3f0)
20.085537f0
```

The exported functions include (within 1 ulp)
```julia
sin, cos, tan, asin, acos, atan, atan2, sincos, sinh, cosh, tanh,
    asinh, acosh, atanh, log, log10, log1p, ilog2, exp, exp2, exp10, expm1, ldexp, cbrt, pow
 ```
 Faster variants include (within 3 ulp)

 ```julia
sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast
```

You can also access `Libm.Musl.log(x)`  for a different implementation of the logarithmic function and `Libm.Musl.erf(x)` and `Libm.Musl.erfc(x)` for the error function and the complementary error function. 

# Benchmarks

You can benchmark the performance of the `Libm.jl` math library on your machine by running
```julia
include(joinpath(Pkg.dir("Libm"), "bench", "bench.jl"))
```
Note this will take some time to run for the first time (subsequent runs will be much faster). Please ensure you run the benchmarks under `-O3` and that you are not concurrently running other cpu intensive tasks.
**Feel free to post your benchmark results at https://github.com/JuliaMath/Libm.jl/issues/34**

Sample benchmark results (commit 1074d479c469e2afcf0d7608ebcbbad7e8aa0a77)
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

sin
time: 0.75 Float64
time: 2.00 Float32

cos
time: 0.82 Float64
time: 2.18 Float32

tan
time: 1.01 Float64
time: 1.94 Float32

asin
time: 2.32 Float64
time: 2.80 Float32

acos
time: 2.96 Float64
time: 3.74 Float32

atan
time: 2.05 Float64
time: 2.30 Float32

exp
time: 1.26 Float64
time: 1.57 Float32

exp2
time: 3.26 Float64
time: 5.69 Float32

exp10
time: 0.25 Float64
time: 0.41 Float32

expm1
time: 2.34 Float64
time: 3.70 Float32

log
time: 1.86 Float64
time: 2.22 Float32

log10
time: 1.61 Float64
time: 1.98 Float32

log1p
time: 1.50 Float64
time: 1.52 Float32

sinh
time: 0.96 Float64
time: 1.02 Float32

cosh
time: 1.30 Float64
time: 1.46 Float32

tanh
time: 1.84 Float64
time: 1.14 Float32

asinh
time: 1.15 Float64
time: 1.27 Float32

acosh
time: 1.09 Float64
time: 1.09 Float32

atanh
time: 2.11 Float64
time: 2.37 Float32

cbrt
time: 2.38 Float64
time: 3.79 Float32

```

