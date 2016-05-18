module SnpArraysTest
using SnpArrays
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end




info("Testing Constructors")
x = rand(0:2, 100, 1000)
snp = SnpArray(x)
@test typeof(snp) <: AbstractMatrix

i, j = rand(1:size(snp, 1)), rand(1:size(snp, 2))

snp[i, j] = 2
@test snp[i, j] == (true, true)

snp[i, j] = NaN
@test snp[i, j] == (false, true)




end # module
