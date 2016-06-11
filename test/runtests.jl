include("../src/SnpArrays.jl")

module SnpArraysTest

using SnpArrays

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end


info("Constructors")
SnpArray("../docs/hapmap3")
gc()
SnpArrayBM("../docs/hapmap3")
gc()
@time sa = SnpArray("../docs/hapmap3")
gc()
@time sabm = SnpArrayBM("../docs/hapmap3")
gc()

info("Indexing")
n,p = size(sa)
for j=1:100
    for i=1:100
        @test sabm[i,j] == sa[i,j]
    end
end

info("summarize")
summarize(sa)
gc()
summarize(sabm)
gc()
@time maf1, _, _, _ = summarize(sa)
gc()
@time maf2, _, _, _ = summarize(sabm)
for i=1:p
    @test maf1[i] == maf2[i]
end
gc()

info("grm")
grm(sa)
gc()
grm(sabm)
gc()
@time grm1 = grm(sa)
gc()
@time grm2 = grm(sabm)
gc()
@show grm1[1,1] grm2[1,1]


info("mom")
grm(sa; method=:MoM)
gc()
grm(sabm; method=:MoM)
gc()
@time mom1 = grm(sa; method=:MoM)
gc()
@time mom2 = grm(sabm; method=:MoM)
gc()
@show mom1[1,1] mom2[1,1]

info("pca")
pca(sa)
gc()
pca(sabm)
gc()
@time _, _, v1 = pca(sa)
gc()
@time _, _, v2 = pca(sabm)
gc()
for i=1:5
    @show v1[i], v2[i]
end


info("pca_sp")
pca_sp(sa)
gc()
pca_sp(sabm)
gc()
@time _, _, v1 = pca_sp(sa)
gc()
@time _, _, v2 = pca_sp(sabm)
gc()
for i=1:5
    @show v1[i], v2[i]
end


#info("grm")
#@time grm(sa; method=:GRM)
#@time grm(sabm; method=:GRM)

#info("mom")
#@time grm(sa; method=:MoM)
#@time grm(sabm; method=:MoM)

#info("pca")
#@time _, _, v1 = pca(sa)

#info("pca_sp")
#@time _, _, v1 = pca_sp(sa)
#@show v1[1]

end # module
