using Libm
using Base.Test

# write your own tests here

x = 2*rand(10000,1)-1
@test Base.erf.(x) ≈ Libm.erf.(x)
@test Base.erfc.(x) ≈ Libm.erfc.(x)