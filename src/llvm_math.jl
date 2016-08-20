using Base: llvmcall

function _llvm_ufun(name, jltype::Type, llvmtype, suffix)
	declarations = """declare $llvmtype @llvm.$name.$suffix($llvmtype)"""
	body = """%2 = call $llvmtype @llvm.$name.$suffix($llvmtype %0)
			   ret $llvmtype %2"""
	quote
		function $name(x::$jltype)
			llvmcall(($declarations, $body),
					$jltype, Tuple{$jltype}, x)
	   end
	end
end


function _llvm_bfun(name, jltype::Type, llvmtype, suffix)
	declarations = """declare $llvmtype @llvm.$name.$suffix($llvmtype, $llvmtype)"""
	body = """%3 = call $llvmtype @llvm.$name.$suffix($llvmtype %0, $llvmtype %1)
			   ret $llvmtype %3"""
	quote
		function $name(x::$jltype, y::$jltype)
			llvmcall(($declarations, $body),
					$jltype, Tuple{$jltype,$jltype}, x, y)
	   end
	end
end



# uniry functions
for fun in [
		:sin,:cos, :exp,:exp2,:log, :log10, :log2,
		:sqrt, :fabs, :ceil, :floor, :trunc, :round, :rint, :nearbyint
		]
	eval(_llvm_ufun(fun, Float64, :double, :f64))
	eval(_llvm_ufun(fun, Float32, :float, :f32))
end


# binary functions
for fun in [:pow, :maxnum, :minnum]
	eval(_llvm_bfun(fun, Float64, :double, :f64))
	eval(_llvm_bfun(fun, Float32, :float, :f32))
end

const fmax=maxnum
const fmin=minnum


#const remainder=rem=frem



