<div align="center"><img src="https://cloud.githubusercontent.com/assets/4319522/19087675/e8a6ea96-8a40-11e6-8f2e-450bef7c1ec9.png" alt="Amal Logo" width="550"></img> </div>

Amal, a pure Julia math library

...with native support for `Float16`, `Float32`, `Float64` types. In the future we plan support for `Double` arithmetic types and `Float80`, `Float128` (when they become available in Julia) types.


The exported functions currently include
```julia
exp, exp2, tan, atan
```
More function to some in the near future. 

*Currently the Amal library is a work in progress*


The Amal library principles include: avoid expensive branches, no table look ups, and also division operations.
We try to adhere to these principles throughout the library in addition to meeting strict accuracy constrains (max 1.5 ulp, but try to be less than < 1 ulp for most functions). To achieve these goals Amal borrows ideas from several open source math libraries, including SLEEF, Cephes, and Musl. As this library develops these guiding principles may adapt to meet different needs.

[![Travis Build Status](https://travis-ci.org/JuliaMath/Amal.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Amal.jl)
[![Appveyor Build Status](https://ci.appveyor.com/api/projects/status/307l6b799amrpvks/branch/master?svg=true)](https://ci.appveyor.com/project/musm/Amal-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaMath/Amal.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/JuliaMath/Amal.jl?branch=master)
[![codecov.io](http://codecov.io/github/JuliaMath/Amal.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaMath/Amal.jl?branch=master)

## Installation

We recommend running julia with `-O3` for maximal performance using `Amal.jl` and to also build a custom system image by running
```julia
julia> is_windows() && (Pkg.add("WinRPM"); using WinRPM; WinRPM.install("gcc"))
julia> include(joinpath(dirname(JULIA_HOME),"share","julia","build_sysimg.jl"))
julia> build_sysimg(force=true)
```
and then to restart `julia`; this will ensure you are taking full advantage of hardware [FMA](https://en.wikipedia.org/wiki/FMA_instruction_set)  if your CPU supports it.

## Usage

To use  `Amal.jl`
```julia
julia> Pkg.clone("https://github.com/JuliaMath/Amal.jl.git")

julia> using Amal

julia> Amal.exp(2.0)
7.38905609893065

julia> Amal.atan(1.0)
0.7853981633974483
```
