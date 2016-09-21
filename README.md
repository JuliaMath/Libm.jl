# Libm.jl

[![Travis Build Status](https://travis-ci.org/JuliaMath/Libm.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Libm.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/307l6b799amrpvks/branch/master?svg=true)](https://ci.appveyor.com/project/simonbyrne/libm-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Libm.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Libm.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Libm.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Libm.jl?branch=master)

This aims to be an implementation of the functions provided by the C math library (aka libm) in pure Julia.

# Usage

We recommend running julia with `-O3` for maximal performance using `Libm.jl` and also to build a custom system image by running
```julia
# Pkg.add("WinRPM"); WinRPM.install("gcc")  # on Windows please first run this line
julia> include(joinpath(dirname(JULIA_HOME),"share","julia","build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then restart your `julia` session. This will ensure you are taking advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set) instruction  if your CPU supports it.


To use  `Libm.jl`
```julia
julia> Pkg.checkout("https://github.com/JuliaMath/Libm.jl.git")
```

Right now most of the useful math functions are accessed with the prefix `x`, for example:
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

The exported functions include
```julia
xatan2, xasin, xacos, xatan, xsin, xcos, xsincos, xtan, xpow, xsinh, xcosh, xtanh,
      xasinh, xacosh, xatanh, xcbrt, xlog, xexp, xexp2, xexp10, xexpm1, xlog10, xlog1p, xilogb, xldexp
 
 # Higher accuracy variants (within 1 ulp). Those not listed below are already within 1 ulp.
 xatan2_u1, xasin_u1, xacos_u1, xatan_u1, xsin_u1, xcos_u1, xsincos_u1, xtan_u1, xcbrt_u1, xlog_u1
```


You can also access `Libm.Musl.log(x)`  for a different implementation of the logarithmic function and `Libm.Musl.erf(x)` and `Libm.Musl.erfc(x)` for the error function and the complementary error function. 