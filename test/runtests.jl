module SnpArraysTest
using SnpArrays
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end


n, p = 100, 1000

info("Constructors")
x = rand(0:2, n, p)
snp = SnpArray(x)
@test typeof(snp) <: AbstractMatrix
@test size(snp) == size(x)
@test eltype(snp) == Tuple{Bool, Bool}


info("getindex and setindex!")
i, j = rand(1:size(snp, 1)), rand(1:size(snp, 2))
snp[i, j] = 2
@test snp[i, j] == (true, true)
snp[i, j] = 1
@test snp[i, j] == (true, false)
snp[i, j] = 0
@test snp[i, j] == (false, false)
snp[i, j] = NaN
@test snp[i, j] == (false, true)


info("convert")
snp_float = Matrix{Float64}(snp)
@test typeof(snp_float) == Matrix{Float64}
@test size(snp) == size(snp_float)
storage = zeros(n, p)
copy!(storage, snp)



end # module
