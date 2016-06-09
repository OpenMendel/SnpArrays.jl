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
#@time sa = SnpArray("/home/huwenbo/Huwenbo/CG10kNhwHg19Clean_v2_Mar2013")
@time sabm = SnpArrayBM("/home/huwenbo/Huwenbo/CG10kNhwHg19Clean_v2_Mar2013")

#info("getindex")
#for i=1:100
#    for j=1:100
#        @test sa[i,j] == sabm[i,j]
#    end
#end

#for i=1:100
#    for j=1:100
#        sa[i,j] = (false,true)
#        @test sa[i,j] == (false,true)
#        sabm[i,j] = (false,true)
#        @test sabm[i,j] == (false,true)
#    end
#end


#info("summarize")
#@time maf, _, _, _ = summarize(sa)
#@time mafbm, _, _, _ = summarize(sabm)
#_, m, n = size(sabm)
#for i=1:n
#    @test maf[i] == mafbm[i]
#end

info("grm")
#@time grm(sa)
@time _grmbm(sabm)

#snp[i, j] = 2
#@test snp[i, j] == (true, true)
#snp[i, j] = 1
#@test snp[i, j] == (true, false)
#snp[i, j] = 0
#@test snp[i, j] == (false, false)
#snp[i, j] = NaN
#@test snp[i, j] == (false, true)


end # module
