# port form musl library, commit e216951f509b71da193da2fc63e25b998740d58b
# see http://git.musl-libc.org/cgit/musl/tree/src/math/erf.c for implementation details
# origin: FreeBSD /usr/src/lib/msun/src/s_erf.c Copyright (c) 1993: Sun Microsystems, Inc.

let
global erf
global erfc

const erx  = 8.45062911510467529297e-01 # 0x3FEB0AC1, 0x60000000

# Coefficients for approximation to  erf on [0,0.84375]
const efx8 =  1.02703333676410069053e+00 # 0x3FF06EBA, 0x8214DB69
const pp0  =  1.28379167095512558561e-01 # 0x3FC06EBA, 0x8214DB68
const pp1  = -3.25042107247001499370e-01 # 0xBFD4CD7D, 0x691CB913
const pp2  = -2.84817495755985104766e-02 # 0xBF9D2A51, 0xDBD7194F
const pp3  = -5.77027029648944159157e-03 # 0xBF77A291, 0x236668E4
const pp4  = -2.37630166566501626084e-05 # 0xBEF8EAD6, 0x120016AC
const qq1  =  3.97917223959155352819e-01 # 0x3FD97779, 0xCDDADC09
const qq2  =  6.50222499887672944485e-02 # 0x3FB0A54C, 0x5536CEBA
const qq3  =  5.08130628187576562776e-03 # 0x3F74D022, 0xC4D36B0F
const qq4  =  1.32494738004321644526e-04 # 0x3F215DC9, 0x221C1A10
const qq5  = -3.96022827877536812320e-06 # 0xBED09C43, 0x42A26120

# Coefficients for approximation to  erf  in [0.84375,1.25]
const pa0  = -2.36211856075265944077e-03 # 0xBF6359B8, 0xBEF77538
const pa1  =  4.14856118683748331666e-01 # 0x3FDA8D00, 0xAD92B34D
const pa2  = -3.72207876035701323847e-01 # 0xBFD7D240, 0xFBB8C3F1
const pa3  =  3.18346619901161753674e-01 # 0x3FD45FCA, 0x805120E4
const pa4  = -1.10894694282396677476e-01 # 0xBFBC6398, 0x3D3E28EC
const pa5  =  3.54783043256182359371e-02 # 0x3FA22A36, 0x599795EB
const pa6  = -2.16637559486879084300e-03 # 0xBF61BF38, 0x0A96073F
const qa1  =  1.06420880400844228286e-01 # 0x3FBB3E66, 0x18EEE323
const qa2  =  5.40397917702171048937e-01 # 0x3FE14AF0, 0x92EB6F33
const qa3  =  7.18286544141962662868e-02 # 0x3FB2635C, 0xD99FE9A7
const qa4  =  1.26171219808761642112e-01 # 0x3FC02660, 0xE763351F
const qa5  =  1.36370839120290507362e-02 # 0x3F8BEDC2, 0x6B51DD1C
const qa6  =  1.19844998467991074170e-02 # 0x3F888B54, 0x5735151D

# Coefficients for approximation to  erfc in [1.25,1/0.35]
const ra0  = -9.86494403484714822705e-03 # 0xBF843412, 0x600D6435 
const ra1  = -6.93858572707181764372e-01 # 0xBFE63416, 0xE4BA7360 
const ra2  = -1.05586262253232909814e+01 # 0xC0251E04, 0x41B0E726 
const ra3  = -6.23753324503260060396e+01 # 0xC04F300A, 0xE4CBA38D 
const ra4  = -1.62396669462573470355e+02 # 0xC0644CB1, 0x84282266 
const ra5  = -1.84605092906711035994e+02 # 0xC067135C, 0xEBCCABB2 
const ra6  = -8.12874355063065934246e+01 # 0xC0545265, 0x57E4D2F2 
const ra7  = -9.81432934416914548592e+00 # 0xC023A0EF, 0xC69AC25C 
const sa1  =  1.96512716674392571292e+01 # 0x4033A6B9, 0xBD707687 
const sa2  =  1.37657754143519042600e+02 # 0x4061350C, 0x526AE721 
const sa3  =  4.34565877475229228821e+02 # 0x407B290D, 0xD58A1A71 
const sa4  =  6.45387271733267880336e+02 # 0x40842B19, 0x21EC2868 
const sa5  =  4.29008140027567833386e+02 # 0x407AD021, 0x57700314 
const sa6  =  1.08635005541779435134e+02 # 0x405B28A3, 0xEE48AE2C 
const sa7  =  6.57024977031928170135e+00 # 0x401A47EF, 0x8E484A93 
const sa8  = -6.04244152148580987438e-02 # 0xBFAEEFF2, 0xEE749A62 

# Coefficients for approximation to  erfc in [1/.35,28]
const rb0  = -9.86494292470009928597e-03 # 0xBF843412, 0x39E86F4A
const rb1  = -7.99283237680523006574e-01 # 0xBFE993BA, 0x70C285DE
const rb2  = -1.77579549177547519889e+01 # 0xC031C209, 0x555F995A
const rb3  = -1.60636384855821916062e+02 # 0xC064145D, 0x43C5ED98
const rb4  = -6.37566443368389627722e+02 # 0xC083EC88, 0x1375F228
const rb5  = -1.02509513161107724954e+03 # 0xC0900461, 0x6A2E5992
const rb6  = -4.83519191608651397019e+02 # 0xC07E384E, 0x9BDC383F
const sb1  =  3.03380607434824582924e+01 # 0x403E568B, 0x261D5190
const sb2  =  3.25792512996573918826e+02 # 0x40745CAE, 0x221B9F0A
const sb3  =  1.53672958608443695994e+03 # 0x409802EB, 0x189D5118
const sb4  =  3.19985821950859553908e+03 # 0x40A8FFB7, 0x688C246A
const sb5  =  2.55305040643316442583e+03 # 0x40A3F219, 0xCEDF3BE6
const sb6  =  4.74528541206955367215e+02 # 0x407DA874, 0xE79FE763
const sb7  = -2.24409524465858183362e+01 # 0xC03670E2, 0x42712D62
    
function erfc1{T<:FloatTypes}(x::T)
    s = abs(x) - 1
    P = @horner_oftype s pa0 pa1 pa2 pa3 pa4 pa5 pa6
    Q = @horner_oftype s 1.0 qa1 qa2 qa3 qa4 qa5 qa6
    return 1 - T(erx) - P/Q
end

function erfc2{T<:FloatTypes}(ix::UInt32, x::T)
    if ix < highword(T(1.25))
        # 0.84375 <= |x| < 1.25
        return erfc1(x)
    end
    # 1.25 <= |x| < 28
    x = abs(x)
    s = 1/(x*x)
    if ix < highword(T(1/0.35000001))
        # 1.25 <= |x| < 1/.35 ~ 2.85714
        R = @horner_oftype s ra0 ra1 ra2 ra3 ra4 ra5 ra6 ra7
        S = @horner_oftype s 1.0 sa1 sa2 sa3 sa4 sa5 sa6 sa7 sa8
    else
        # 1/.35 <= |x| < 28
        R = @horner_oftype s rb0 rb1 rb2 rb3 rb4 rb5 rb6
        S = @horner_oftype s 1.0 sb1 sb2 sb3 sb4 sb5 sb6 sb7
    end
    z = trunclo(x)
    return exp(-z*z-T(0.5625))*exp((z-x)*(z+x)+R/S)/x
end

function erf{T<:FloatTypes}(x::T)
    ix = highword(x)
    sign = (ix>>31) % Bool
    ix &= 0x7fffffff
    if ix >= highword(T(Inf)) # erf(nan)=nan, erf(+-inf)=+-1
        return 1-2*sign + 1/x
    end
    if ix < highword(T(0.84375)) # |x| < 0.84375
        if ix < highword(T(2)^-28) #|x| < 2**-28  avoid underflow
            return (8*x +T(efx8)*x)/8
        end
        z = x*x
        r = @horner_oftype z pp0 pp1 pp2 pp3 pp4
        s = @horner_oftype z 1 qq1 qq2 qq3 qq4 qq5
        y = r/s
        return x + x*y
    end
    if ix < highword(T(6)) # 0.84375 <= |x| < 6
        y = 1 - erfc2(ix,x)
    else
        y = 1 - realmin(T)
    end
    return sign ? -y : y
end

function erfc{T<:FloatTypes}(x::T)
    ix = highword(x)
    sign = (ix>>31) % Bool
    ix &= 0x7fffffff
    if ix >= highword(T(Inf)) # erfc(nan)=nan, erfc(+-inf)=0,2
        return 2*sign + 1/x
    end
    if ix < highword(T(0.84375)) # |x| < 0.84375
        if ix < highword(T(2)^-56)  # |x| < 2**-56
            return 1 - x
        end
        z = x*x
        r = @horner_oftype z pp0 pp1 pp2 pp3 pp4
        s = @horner_oftype z 1 qq1 qq2 qq3 qq4 qq5
        y = r/s
        if sign || ix < highword(T(0.25)) # x < 1/4 
            return 1 - (x+x*y)
        end
        return T(0.5) - (x - T(0.5) + x*y)
    end
    if ix < highword(T(28)) # 0.84375 <= |x| < 28
        return sign ? 2 - erfc2(ix,x) : erfc2(ix,x)
    end
    return sign ? 2 - realmin(T) : realmin(T)*realmin(T)
end

end
