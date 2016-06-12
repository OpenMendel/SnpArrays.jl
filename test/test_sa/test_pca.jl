include("../../src/SnpArrays.jl")

module SnpArraysTest

using SnpArrays

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

sa = SnpArray("/home/huwenbo/Huwenbo/merge-geno")
gc()
@time pca(sa)
gc()

end
