using Libm
using Suppressor
using Base.Test

using Base: significand_bits
function countulp(T,X::AbstractFloat,Y::AbstractFloat)
    x,y = T(X),T(Y)
    (isnan(x) && isnan(y)) && return 0
    (isnan(x) || isnan(y)) && return 10000
    (isinf(x) && isinf(y)) && return (signbit(x) == signbit(y)) ? 0 : 10001
    (x == 0 && y == 0) && return (signbit(x) == signbit(y)) ? 0 : 10002
    if isfinite(x) && isfinite(y)
        d = abs(X - Y)
        return T(d/eps(y))
        # k = abs(frexp(Y)[2])
        # return T(ldexp(d,significand_bits(T)-k+1))
    end
    return 10003
end
function cmpdenorm{Tx<:AbstractFloat, Ty<:AbstractFloat}(x::Tx, y::Ty)
    (isnan(x) && isnan(y)) && return signbit(x) == signbit(y)
    (isinf(x) && isinf(y)) && return signbit(x) == signbit(y)
    (isfinite(x) && isfinite(y)) && return signbit(x) == signbit(y)
    return false
end



# # the following compares the ulp between x and y.
# # First it promotes them to the larger of the two types x,y
# infh(::Type{Float64}) = 1e+300
# infh(::Type{Float32}) = 1e+37
# function countulp(T, x::AbstractFloat, y::AbstractFloat)
#     X,Y = promote(x,y)
#     x, y = T(X), T(Y) # Cast to smaller type
#     (isnan(x) && isnan(y)) && return 0
#     (isnan(x) || isnan(y)) && return 10000
#     if isinf(x)
#         (sign(x) == sign(y) && abs(y) > infh(T)) && return 0 # Relaxed infinity handling
#         return 10001
#     end
#     (x ==  Inf && y ==  Inf) && return 0
#     (x == -Inf && y == -Inf) && return 0
#     if y == 0
#         (x == 0) && return 0
#         return 10002
#     end
#     if isfinite(x) && isfinite(y)
#         return T(abs(X - Y)/eps(y))
#     end
#     return 10003
# end
# countulp{T<:AbstractFloat}(x::T, y::T) = countulp(T,x,y)


# overide domain checking that base adheres to
using Base.MPFR.ROUNDING_MODE
for f in (:sin, :cos, :tan, :asin, :acos, :atan, :asinh, :acosh, :atanh, :log, :log10, :log2, :log1p)
    @suppress @eval begin
        import Base.$f
        function ($f)(x::BigFloat)
            z = BigFloat()
            ccall($(string(:mpfr_,f), :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32), &z, &x, ROUNDING_MODE[])
            return z
        end
    end
end

strip_module_name(f::Function) = last(split(string(f), '.')) # strip module name from function f

# test the accuracy of a function where fun_table is a Dict mapping the function you want
# to test to a reference function
# xx is an array of values (which may be tuples for multiple arugment functions)
# tol is the acceptable tolerance to test against
function test_acc(T, fun_table, xx, tol; debug=false, tol_debug=5)
    @testset "accuracy $(strip_module_name(xfun))" for (xfun, fun) in fun_table
        rmax = 0.0
        rmean = 0.0
        xmax = map(zero, first(xx))
        for x in xx
            q = xfun(x...)
            c = fun(map(BigFloat,x)...)
            u = countulp(T,q,c)
            rmax = max(rmax, u)
            xmax = rmax == u ? x : xmax
            rmean += u
            if debug && u > tol_debug
                @printf("%s = %.20g\n%s  = %.20g\nx = %.20g\nulp = %g\n", strip_module_name(xfun), q, strip_module_name(fun), T(c), x, ulp(T(c)))
            end
        end
        rmean = rmean/length(xx)

        t = @test trunc(rmax,1) <= tol

        fmtxloc = isa(xmax, Tuple) ? string('(', join((@sprintf("%.5g", x) for x in xmax), ", "), ')') : @sprintf("%.5f", xmax)
        println(rpad(strip_module_name(xfun), 18, " "), ": max ", @sprintf("%g", rmax),
            rpad(" at x = "*fmtxloc, 40, " "),
            ": mean ", @sprintf("%g", rmean))
    end
end

function runtests()
    include("dnml_nan.jl")
    include("accuracy.jl")
    include("accuracy_test.jl")
end

runtests()
