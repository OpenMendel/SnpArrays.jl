using SnpArrays, BenchmarkTools

hapmap = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3")
n, snps = size(hapmap)

outv = zeros(n)
# @code_warntype copy!(outv, view(hapmap, :, 1), :additive, false, false, false)
# @benchmark copy!(outv, view(hapmap, :, 1), :additive, false, false, false)
# @benchmark copy!(outv, hapmap[:, 1], :additive, false, false, false)
@time for s in 1:snps
    copy!(outv, view(hapmap, :, s), :additive, false, false, false)
end

# outm = zeros(n, snps)
# @time copy!(outm, hapmap)
# @benchmark copy!(outm, hapmap)
# copy!(outm, hapmap, :additive, false, false, false)
# @code_warntype copy!(outm, hapmap, :additive, false, false, false)

# a = (false, true)
# @code_warntype isnan(a)
# @code_llvm isnan(a)

# SnpArrays.convert(Float64, a, true, :additive)
# @code_warntype SnpArrays.convert(Float64, a, true, :additive)
# @benchmark SnpArrays.convert(Float64, a, true, :additive)
