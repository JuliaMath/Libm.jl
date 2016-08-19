using Libm
using Base.Test

# write your own tests here

x = 2*rand(10000,1)-1
@test Base.erf.(x) ≈ Libm.erf.(x)
@test isnan(Base.erf.(NaN)) 
@test Base.erf.(Inf) == Libm.erf.(Inf)
@test Base.erf.(-Inf) == Libm.erf.(-Inf)

@test Base.erfc.(x) ≈ Libm.erfc.(x)
@test isnan(Base.erfc.(NaN))
@test Base.erfc.(Inf) == Libm.erfc.(Inf)
@test Base.erfc.(-Inf) == Libm.erfc.(-Inf)