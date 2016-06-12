include("../../src/SnpArrays.jl")

module SnpArraysTest

using SnpArrays

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

sa = SnpArrayBM("/home/huwenbo/Huwenbo/merge-geno")
_,n,p = size(sa)
v = vec(rand(p,1))
gc()
@time sa*v
gc()

end
