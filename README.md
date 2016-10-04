<div align="center"><img src="https://cloud.githubusercontent.com/assets/4319522/19087675/e8a6ea96-8a40-11e6-8f2e-450bef7c1ec9.png" alt="Amal Logo" width="550"></img> </div>

Amal, a pure Julia math library

...with native support for `Float16`, `Float32`, `Float64` types. In the future we plan support for `Double` arithmetic types and `Float80`, `Float128` types.


The Amal library principles include: avoid expensive branches, no table look ups, and avoid expensive division operations.
We try to adhere to these principles throughout the library in addition to meeting strict accuracy constrains (max 1.5 ulp, but try to be less than < 1 ulp for most functions). To achieve these goals Amal borrows ideas from several open source math libraries.

[![Travis Build Status](https://travis-ci.org/JuliaMath/Libm.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Libm.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/307l6b799amrpvks/branch/master?svg=true)](https://ci.appveyor.com/project/simonbyrne/libm-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Libm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Libm.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Libm.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Libm.jl?branch=master)

## Installation

We recommend running julia with `-O3` for maximal performance using `Libm.jl` and to also build a custom system image by running
```julia
# Pkg.add("WinRPM"); WinRPM.install("gcc")  # on Windows please first run this line
julia> include(joinpath(dirname(JULIA_HOME),"share","julia","build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then to restart `julia`; this will ensure you are taking full advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set)  if your CPU supports it.

## Usage

To use  `Libm.jl`
```julia
julia> Pkg.clone("https://github.com/JuliaMath/Libm.jl.git")
```

Currently, the following functions are available in the Amal library: `exp, exp2, atan`

```julia
julia> using Libm.Amal

julia> Amal.exp(2.0)
7.38905609893065

julia> Amal.atan(1.0)
0.7853981633974483
```


## The Sleef submodule

The `Sleef` submodule is a full feature complete port (including additional micro optimzations) of the sleef C library.

The `Sleef` submodule functions include (within 1 ulp)
```julia
sin, cos, tan, asin, acos, atan, atan2, sincos, sinh, cosh, tanh,
    asinh, acosh, atanh, log, log2, log10, log1p, ilog2, exp, exp2, exp10, expm1, ldexp, cbrt, pow
 ```
 Faster variants include (within 3 ulp)

 ```julia
sin_fast, cos_fast, tan_fast, sincos_fast, asin_fast, acos_fast, atan_fast, atan2_fast, log_fast, cbrt_fast
```

These function can be accessed as follows
```julia
julia> using Libm.Sleef

julia> Sleef.exp(2.0)
7.38905609893065
```


### Benchmarks

<!-- You can benchmark the performance of the `Libm.jl` math library on your machine by running
```julia
include(joinpath(Pkg.dir("Libm"), "benchmark", "benchmark.jl"))
```
Note this will take some time to run for the first time (subsequent runs will be much faster). Please ensure you run the benchmarks under `-O3` and that you are not concurrently running other cpu intensive tasks.
**Feel free to post your benchmark results at https://github.com/JuliaMath/Libm.jl/issues/34**
 -->

Sample benchmark results (commit 053528ab9b9e673bca99d1dc7db5445e39ad0d77) of the Sleef submodule

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

Benchmarks: median ratio Libm/Base
sin                
time: 0.63 Float64 
time: 1.69 Float32 
                   
cos                
time: 0.71 Float64 
time: 1.73 Float32 
                   
tan                
time: 0.44 Float64 
time: 1.53 Float32 
                   
asin               
time: 1.82 Float64 
time: 2.38 Float32 
                   
acos               
time: 2.82 Float64 
time: 2.66 Float32 
                   
atan               
time: 1.59 Float64 
time: 1.91 Float32 
                   
exp                
time: 0.67 Float64 
time: 0.96 Float32 
                   
exp2               
time: 1.77 Float64 
time: 4.80 Float32 
                   
exp10              
time: 0.21 Float64 
time: 0.30 Float32 
                   
expm1              
time: 1.95 Float64 
time: 2.27 Float32 
                   
log                
time: 1.66 Float64 
time: 2.01 Float32 
                   
log2               
time: 1.37 Float64 
time: 2.04 Float32 
                   
log10              
time: 0.82 Float64 
time: 1.22 Float32 
                   
log1p              
time: 1.07 Float64 
time: 1.21 Float32 
                   
sinh               
time: 0.67 Float64 
time: 0.86 Float32 
                   
cosh               
time: 0.78 Float64 
time: 0.85 Float32 
                   
tanh               
time: 0.93 Float64 
time: 0.98 Float32 
                   
asinh              
time: 1.08 Float64 
time: 1.06 Float32 
                   
acosh              
time: 0.62 Float64 
time: 0.95 Float32 
                   
atanh              
time: 1.54 Float64 
time: 1.92 Float32 
                   
cbrt               
time: 2.49 Float64 
time: 3.75 Float32 

```

