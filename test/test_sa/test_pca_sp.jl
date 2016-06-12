include("../../src/SnpArrays.jl")

module SnpArraysTest

using SnpArrays

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

sa = SnpArray("../../docs/hapmap3")
#sa = SnpArray("/home/huwenbo/Huwenbo/merge-geno")
#sa = SnpArray(sa.A1[:,1:25:end], sa.A2[:,1:25:end])
gc()
@time _,_,v = pca_sp(sa)
gc()
@show v[1]

end
