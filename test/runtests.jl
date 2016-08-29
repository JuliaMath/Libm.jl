using Libm
using Base.Test

isnzero{T<:AbstractFloat}(x::T) = signbit(x)
ispzero{T<:AbstractFloat}(x::T) = !signbit(x)

function cmpdenorm{Tx<:AbstractFloat, Ty<:AbstractFloat}(x::Tx, y::Ty)
    sizeof(Tx) < sizeof(Ty) ? y = Tx(y) : x = Ty(x) # cast larger type to smaller type
    (isnan(x) && isnan(y)) && return true
    (isnan(x) || isnan(y)) && return false
    (isinf(x) != isinf(y)) && return false
    (x == Tx(Inf)  && y == Ty(Inf))  && return true
    (x == Tx(-Inf) && y == Ty(-Inf)) && return true
    if y == 0
        (ispzero(x) && ispzero(y)) && return true
        (isnzero(x) && isnzero(y)) && return true
        return false
    end
    (!isnan(x) && !isnan(y) && !isinf(x) && !isinf(y)) && return sign(x) == sign(y)
    return false
end

# the following compares x with y first it promotes them to the larger of the two types x,y
# todo update eps to ulp from reference article
countulp{T<:Float64}(x::T, y::BigFloat) = _countulp(T, promote(x, y)...)
countulp{T<:Float32}(x::T, y::Union{Float64, BigFloat})  = _countulp(T, promote(x, y)...)
function _countulp(T, x::AbstractFloat, y::AbstractFloat)::T
    T == Float64 && (infh = 1e+300)
    T == Float32 && (infh = 1e+37)
    fx, fy = T(x), T(y) # cast to smaller type
    (isnan(fx) && isnan(fy)) && return 0
    (isnan(fx) || isnan(fy)) && return 10000
    if isinf(fx)
        if sign(fx) == sign(fy) && abs(fy) > infh
            return 0 #Relaxed infinity handling
        end
        return 10001
    end
    (fx ==  Inf && fy ==  Inf) && return 0
    (fx == -Inf && fy == -Inf) && return 0
    if fy == 0
        if fx == 0
            return 0
        end
    return 10002
    end
    if !isnan(fx) && !isnan(fy) && !isinf(fx) && !isinf(fy)
        return T(abs( (x - y) / eps(T(y)) ))
    end
    return 10003
end

# get rid off annoying warnings from overwritten function
macro nowarn(expr)
    quote
        stderr = STDERR
        stream = open("null", "w")
        redirect_stderr(stream)
        result = $(esc(expr))
        redirect_stderr(stderr)
        close(stream)
        result
    end
end

# overide domain checking that base adheres to
using Base.MPFR.ROUNDING_MODE
for f in (:sin, :cos, :tan, :asin, :acos, :atan, :asinh, :acosh, :atanh, :log, :log10, :log1p)
    @eval begin
        import Base.$f
        @nowarn function ($f)(x::BigFloat)
            z = BigFloat()
            ccall($(string(:mpfr_,f), :libmpfr), Int32, (Ptr{BigFloat}, Ptr{BigFloat}, Int32), &z, &x, ROUNDING_MODE[])
            return z
        end
    end
end

strip_module_name(f::Function) = split(string(f), '.')[end] # strip module name from function f
"""
test the accuracy of a function where fun_table is a Dict mapping the function you want
to test to a reference function
xx is an array of values (which may be tuples for multiple arugment functions)
tol is the acceptable tolerance to test against
"""
function test_acc(T, fun_table, xx, tol; debug=false, tol_debug=5)
    @testset "accuracy $(strip_module_name(xfun))" for (xfun, fun) in fun_table
        rmax = 0.0
        for x in xx
            q = xfun(x...)
            c = fun(map(BigFloat,x)...)
            u = countulp(q, c)
            rmax = max(rmax, u)
            if debug && rmax > tol_debug
                @printf("%s = %.20g\n%s  = %.20g\nx = %.20g\nulp = %g\n", strip_module_name(xfun), q, strip_module_name(fun), T(c), x, ulp(T,c))
            end
        end
        t = @test rmax < tol
        t.value == true ? (v = "GOOD") : (v = "FAIL")
        println(rpad(strip_module_name(xfun), 15, " "), " : ", @sprintf("%f", rmax), " ... ", v)
    end
end

const pow = ^
function runtests()
    include("dnml_nan.jl")
    include("accuracy.jl")
    # include("accuracy_base.jl") # uncomment to benchmark base
end

runtests()

