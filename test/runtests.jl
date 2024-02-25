using LinearAlgebra, SnpArrays, SparseArrays, Test

const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed")) # no missing genotypes
const mouse = SnpArray(SnpArrays.datadir("mouse.bed")) # has missing genotypes

@testset "size" begin
@test size(EUR) == (379, 54051)
@test size(EUR, 1) == 379
@test size(EUR, 2) == 54051
@test size(mouse) == (1940, 10150)
@test size(mouse, 3) == 1
@test_throws ErrorException size(mouse, 0)
@test length(mouse) == 1940 * 10150
@test eltype(mouse) == UInt8
end

@testset "counts" begin
cc = counts(EUR, dims=1)
@test counts(EUR) == sum(cc, dims=2)
@test size(cc) == (4, size(EUR, 2))
@test view(cc, :, 1) == [2, 0, 70, 307]
@test all(sum(cc, dims = 1) .== size(EUR, 1))
@test all(iszero, view(cc, 2, :))
rc = counts(EUR, dims=2)
@test size(rc) == (4, size(EUR, 1))
@test view(rc, :, 1) == [2997, 0, 13143, 37911]
@test all(sum(rc, dims = 1) .== size(EUR, 2))
@test all(iszero, view(rc, 2, :))
@test_throws ArgumentError counts(EUR, dims=3)
end

@testset "means" begin
# ADDITIVE_MODEL
cmns = mean(EUR, dims = 1)
@test size(cmns) == (1, size(EUR, 2))
@test cmns[1] ≈ 1.804749340369393
@test mean(cmns) ≈ mean(EUR)
@test minimum(cmns) ≈ 1.0
@test maximum(cmns) ≈ 1.9788918205804749
rmns = mean(EUR, dims = 2)
@test size(rmns) == (size(EUR, 1), 1)
@test rmns[1] ≈ 1.6459454959205195
@test mean(rmns) ≈ mean(EUR)
rmnextrma = extrema(rmns)
@test rmnextrma[1] ≈ 1.637749532848606
@test rmnextrma[2] ≈ 1.6627259440158368
@test_throws ArgumentError mean(EUR, dims=3)
# DOMINANT_MODEL
cmns = mean(EUR, dims = 1, model = DOMINANT_MODEL)
@test size(cmns) == (1, size(EUR, 2))
@test cmns[1] ≈ 0.9947229551451188
@test mean(cmns) ≈ mean(EUR, model = DOMINANT_MODEL)
@test minimum(cmns) ≈ 0.5989445910290238
@test maximum(cmns) ≈ 1.0
rmns = mean(EUR, dims = 2, model = DOMINANT_MODEL)
@test size(rmns) == (size(EUR, 1), 1)
@test rmns[1] ≈ 0.9445523672087472
@test mean(rmns) ≈ mean(EUR, model = DOMINANT_MODEL)
rmnextrma = extrema(rmns)
@test rmnextrma[1] ≈ 0.9325266877578583
@test rmnextrma[2] ≈ 0.9530073449149877
@test_throws ArgumentError mean(EUR, dims=3, model=DOMINANT_MODEL)
# RECESSIVE_MODEL
cmns = mean(EUR, dims = 1, model = RECESSIVE_MODEL)
@test size(cmns) == (1, size(EUR, 2))
@test cmns[1] ≈ 0.8100263852242744
@test mean(cmns) ≈ mean(EUR, model = RECESSIVE_MODEL)
@test minimum(cmns) ≈ 0.0
@test maximum(cmns) ≈ 0.9868073878627969
rmns = mean(EUR, dims = 2, model = RECESSIVE_MODEL)
@test size(rmns) == (size(EUR, 1), 1)
@test rmns[1] ≈ 0.7013931287117722
@test mean(rmns) ≈ mean(EUR, model = RECESSIVE_MODEL)
rmnextrma = extrema(rmns)
@test rmnextrma[1] ≈ 0.6932341677304767
@test rmnextrma[2] ≈ 0.7191171301178516
@test_throws ArgumentError mean(EUR, dims=3, model=RECESSIVE_MODEL)
end

@testset "var/std" begin
cvars = var(EUR, dims = 1)
@test size(cvars) == (1, size(EUR, 2))
@test cvars[1] ≈ 0.16812553224162724
@test iszero(minimum(cvars))
@test maximum(cvars) ≈ 0.8705727966941688
@test std(EUR, dims = 1) ≈ sqrt.(cvars)
rvars = var(EUR, dims = 2)
@test size(rvars) == (size(EUR, 1), 1)
@test rvars[1] ≈ 0.33960146078503206
@test minimum(rvars) ≈ 0.31972707244268267
@test maximum(rvars) ≈ 0.365901927927013
@test std(EUR, dims = 2) ≈ sqrt.(rvars)
end

@testset "missingpos" begin
mp = missingpos(mouse)
cc = counts(mouse, dims=1)
@test isa(mp, SparseMatrixCSC)
@test sum(mp, dims = 1) == view(cc, 2:2, :)
@test sum(mp, dims = 2) == view(counts(mouse, dims=2), 2:2, :)'
end

@testset "create bed" begin
tmpbf = SnpArray("tmp.bed", 5, 3)
@test isfile("tmp.bed")
@test all(tmpbf .== 0x00)
fill!(tmpbf, 0x02)
tmpbf2 = SnpArray("tmp.bed", 5)
@test all(tmpbf2 .== 0x02)
Sys.iswindows() || rm("tmp.bed", force=true)
tmpbf2 = SnpArray("tmp2.bed", SnpArray(undef, 5, 3))
fill!(tmpbf2, 0x01)
@test all(tmpbf2 .== 0x01)
Sys.iswindows() || rm("tmp2.bed", force=true)
end

@testset "convert" begin
# check convert coding
tmpbf = SnpArray(undef, 4, 1)
tmpbf[1] = 0x00
tmpbf[2] = 0x01
tmpbf[3] = 0x02
tmpbf[4] = 0x03
# additive model
v = convert(Matrix{Float64}, tmpbf, model=ADDITIVE_MODEL)
@test all([v[1], v[3], v[4]] .== [0.0, 1.0, 2.0])
@test isnan(v[2])
# dominant model
v = convert(Matrix{Float64}, tmpbf, model=DOMINANT_MODEL)
@test all([v[1], v[3], v[4]] .== [0.0, 1.0, 1.0])
@test isnan(v[2])
# recessive model
v = convert(Matrix{Float64}, tmpbf, model=RECESSIVE_MODEL)
@test all([v[1], v[3], v[4]] .== [0.0, 0.0, 1.0])
@test isnan(v[2])
# mouse test data
cmns = mean(mouse, dims = 1)
cvrs = var(mouse, dims = 1)
cstd = std(mouse, dims = 1)
v = convert(Vector{Float64}, @view(mouse[:, 1]))
@test isnan(mean(v))
@test mean(filter(!isnan, v)) ≈ cmns[1]
@test isnan(var(v))
@test var(filter(!isnan, v)) ≈ cvrs[1]
@test isnan(std(v))
@test std(filter(!isnan, v)) ≈ cstd[1]
end

@testset "convert (ADMIXTURE)" begin
# convert A2 allele frequency to dosage
@test convert(Float64, 0.5, ADDITIVE_MODEL)  == 1.0
@test convert(Float64, 0.5, DOMINANT_MODEL)  == 0.75
@test convert(Float64, 0.5, RECESSIVE_MODEL) == 0.25
end

@testset "maf" begin
cc = counts(mouse, dims=1)
@test all(0 .≤ maf(mouse) .≤ 0.5)
@test all(minorallele(mouse) .== (cc[1, :] .> cc[4, :]))
cc = counts(EUR, dims=1)
@test all(0 .≤ maf(EUR) .≤ 0.5)
@test all(minorallele(EUR) .== (cc[1, :] .> cc[4, :]))
end

@testset "grm" begin
Φgrm = grm(EUR, method=:GRM)
@test size(Φgrm) == (size(EUR, 1), size(EUR, 1))
@test issymmetric(Φgrm)
@test eigmin(Φgrm) > -1e-8
Φmom = grm(EUR, method=:MoM)
@test size(Φmom) == (size(EUR, 1), size(EUR, 1))
@test issymmetric(Φmom)
Φrbs = grm(EUR, method=:Robust)
@test size(Φrbs) == (size(EUR, 1), size(EUR, 1))
@test issymmetric(Φrbs)
@test eigmin(Φrbs) > -1e-8
m, n, K = size(EUR, 1), size(EUR, 2), 3
P = rand(n, K)
Q = rand(m, K)
Q ./= sum(Q, dims=2)
Φreap = grm_admixture(EUR, transpose(P), transpose(Q))
@test size(Φreap) == (size(EUR, 1), size(EUR, 1))
@test issymmetric(Φreap)
@test eigmin(Φrbs) > -1e-8
end

@testset "filter" begin
rowmask, colmask =  SnpArrays.filter(mouse, min_success_rate_per_row=0.99,
    min_success_rate_per_col=0.99)
SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask; des="tmp")
tmpbf = SnpArray("tmp.bed")
@test size(tmpbf) == (1912, 9997)
@test all(missingrate(tmpbf, 1) .≤ 0.01)
@test all(missingrate(tmpbf, 2) .≤ 0.01)
Sys.iswindows() || rm("tmp.bed", force=true)
rm("tmp.bim", force=true)
rm("tmp.fam", force=true)
rowmask, colmask =  SnpArrays.filter(mouse, min_success_rate_per_row=0.99,
    min_success_rate_per_col=0.99, min_maf=0.01, min_hwe_pval=1e-8)
@test (count(rowmask), count(colmask)) == (1911, 9510)
compress_plink(SnpArrays.datadir("mouse"), "gz")
SnpArrays.filter(SnpArrays.datadir("mouse.bed.gz"), SnpArrays.datadir("mouse.bim.gz"),
SnpArrays.datadir("mouse.fam.gz"), 1:5, 1:3; des="tmp")
tmpbf = SnpArray("tmp.bed")
@test size(tmpbf) == (5, 3)
@test isfile("tmp.bed")
@test isfile("tmp.fam")
@test isfile("tmp.bim")
rm(SnpArrays.datadir("mouse.bed.gz"), force=true)
rm(SnpArrays.datadir("mouse.fam.gz"), force=true)
rm(SnpArrays.datadir("mouse.bim.gz"), force=true)
Sys.iswindows() || rm("tmp.bed", force=true)
rm("tmp.bim", force=true)
rm("tmp.fam", force=true)
end

@testset "SnpBitMatrix-vector multiplication" begin
reltol = 5e-4
for t in [Float32, Float64]
    v1 = randn(t, size(EUR, 1))
    v2 = randn(t, size(EUR, 2))
    for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
        @test norm(SnpBitMatrix{t}(EUR, model=model) * v2 -
            convert(Matrix{t}, EUR, model=model) * v2) /
            norm(convert(Matrix{t}, EUR, model=model) * v2) < reltol
        @test norm(SnpBitMatrix{t}(EUR, center=true, model=model) * v2 -
            convert(Matrix{t}, EUR, center=true, model=model) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, model=model) * v2) < reltol
        @test norm(SnpBitMatrix{t}(EUR, scale=true, model=model) * v2 -
            convert(Matrix{t}, EUR, scale=true, model=model) * v2) /
            norm(convert(Matrix{t}, EUR, scale=true, model=model) * v2) < reltol
        @test norm(SnpBitMatrix{t}(EUR, center=true, scale=true, model=model) * v2 -
            convert(Matrix{t}, EUR, center=true, scale=true, model=model) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, scale=true, model=model) * v2) < reltol
        @test norm(transpose(SnpBitMatrix{t}(EUR, model=model)) * v1 -
            transpose(convert(Matrix{t}, EUR, model=model)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, model=model)) * v1) < reltol
        @test norm(transpose(SnpBitMatrix{t}(EUR, center=true, model=model)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, model=model)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, model=model)) * v1) < reltol
        @test norm(transpose(SnpBitMatrix{t}(EUR, scale=true, model=model)) * v1 -
            transpose(convert(Matrix{t}, EUR, scale=true, model=model)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, scale=true, model=model)) * v1) < reltol
        @test norm(transpose(SnpBitMatrix{t}(EUR, center=true, scale=true, model=model)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model)) * v1) < reltol
    end
end
end

@testset "copyto bitmatrix" begin
for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL], t in [Float32, Float64]
    mousebm = SnpBitMatrix{t}(mouse, model=ADDITIVE_MODEL, center=false, scale=false)
    v = copyto!(zeros(1), @view(mousebm[702]))
    @test v[1] == 0
    mousebm = SnpBitMatrix{t}(mouse, model=ADDITIVE_MODEL, center=true, scale=true)
    v = copyto!(zeros(1), @view(mousebm[702]))
    @test v[1] ≈ -1 * mousebm.μ[1] * mousebm.σinv[1]
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=true, scale=true)
    EURtrue = convert(Matrix{t}, EUR, model=model, center=true, scale=true)
    EURtest = copyto!(zeros(t, size(EUR)), EURbm)
    @test all(EURtest .≈ EURtrue)
    EURtest = copyto!(zeros(10, 5), @view(EURbm[11:20, 1:2:10]))
    @test all(EURtest .≈ EURtrue[11:20, 1:2:10])
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=false, scale=false)
    EURtrue = convert(Matrix{t}, EUR, model=model, center=false, scale=false)
    EURtest = copyto!(zeros(t, size(EUR)), EURbm)
    @test all(EURtest .≈ EURtrue)
    EURtest = copyto!(zeros(10, 5), @view(EURbm[11:20, 1:2:10]))
    @test all(EURtest .≈ EURtrue[11:20, 1:2:10])
end
end

@testset "convert (SnpBitMatrix)" begin
for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL], t in [Float32, Float64]
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=true, scale=true)
    xtrue = convert(Matrix{t}, EUR, model=model, center=true, scale=true)
    xtest = convert(Matrix{t}, EURbm)
    @test all(xtrue .≈ xtest)
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=false, scale=false)
    xtrue = convert(Matrix{t}, EUR, model=model, center=false, scale=false)
    xtest = convert(Matrix{t}, EURbm)
    @test all(xtrue .≈ xtest)
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=true, scale=true)
    xtrue = convert(Matrix{t}, @view(EUR[:, 2:2:10]), model=model, center=true, scale=true)
    xtest = convert(Matrix{t}, @view(EURbm[:, 2:2:10]))
    @test all(xtrue .≈ xtest)
    EURbm = SnpBitMatrix{t}(EUR, model=model, center=false, scale=false)
    xtrue = convert(Matrix{t}, @view(EUR[:, 2:2:10]), model=model, center=false, scale=false)
    xtest = convert(Matrix{t}, @view(EURbm[:, 2:2:10]))
    @test all(xtrue .≈ xtest)
end
end

@testset "copyto SnpLinAlg" begin
for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL], t in [Float32, Float64]
    # imputing missing data
    mousela = SnpLinAlg{t}(mouse, model=ADDITIVE_MODEL, center=true, scale=true, impute=true)
    v = copyto!(zeros(1), @view(mousela[702])) # missing entry
    @test isapprox(v[1], 1.113003134727478, atol=1e-6)
    # not imputing missing data
    mousela = SnpLinAlg{t}(mouse, model=ADDITIVE_MODEL, center=true, scale=true, impute=false)
    v = copyto!(zeros(1), @view(mousela[702])) # missing entry
    @test isnan(v[1])
    # no missing data
    EURla = SnpLinAlg{t}(EUR, model=model, center=true, scale=true)
    EURtrue = convert(Matrix{t}, EUR, model=model, center=true, scale=true)
    EURtest = copyto!(zeros(t, size(EUR)), EURla)
    @test all(EURtest .≈ EURtrue)
    EURtest = copyto!(zeros(10, 5), @view(EURla[11:20, 1:2:10]))
    @test all(EURtest .≈ EURtrue[11:20, 1:2:10])
    EURla = SnpLinAlg{t}(EUR, model=model, center=false, scale=false)
    EURtrue = convert(Matrix{t}, EUR, model=model, center=false, scale=false)
    EURtest = copyto!(zeros(t, size(EUR)), EURla)
    @test all(EURtest .≈ EURtrue)
    EURtest = copyto!(zeros(10, 5), @view(EURla[11:20, 1:2:10]))
    @test all(EURtest .≈ EURtrue[11:20, 1:2:10])
end
end

@testset "convert (SnpLinAlg)" begin
for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL], t in [Float32, Float64]
    EURla = SnpLinAlg{t}(EUR, model=model, center=true, scale=true)
    xtrue = convert(Matrix{t}, EUR, model=model, center=true, scale=true)
    xtest = convert(Matrix{t}, EURla)
    @test all(xtrue .≈ xtest)
    EURla = SnpLinAlg{t}(EUR, model=model, center=false, scale=false)
    xtrue = convert(Matrix{t}, EUR, model=model, center=false, scale=false)
    xtest = convert(Matrix{t}, EURla)
    @test all(xtrue .≈ xtest)
    EURla = SnpLinAlg{t}(EUR, model=model, center=true, scale=true)
    xtrue = convert(Matrix{t}, @view(EUR[:, 2:2:10]), model=model, center=true, scale=true)
    xtest = convert(Matrix{t}, @view(EURla[:, 2:2:10]))
    @test all(xtrue .≈ xtest)
    EURla = SnpLinAlg{t}(EUR, model=model, center=false, scale=false)
    xtrue = convert(Matrix{t}, @view(EUR[:, 2:2:10]), model=model, center=false, scale=false)
    xtest = convert(Matrix{t}, @view(EURla[:, 2:2:10]))
    @test all(xtrue .≈ xtest)
end
end

@testset "SnpLinAlg-vector multiplication (zeroimpute)" begin
reltol = 5e-4
for t in [Float32, Float64]
    v1 = randn(t, size(EUR, 1))
    v2 = randn(t, size(EUR, 2))
    for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
        @test norm(SnpLinAlg{t}(EUR, model=model, impute=false) * v2 -
            convert(Matrix{t}, EUR, model=model, impute=false) * v2) /
            norm(convert(Matrix{t}, EUR, model=model, impute=false) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, center=true, model=model, impute=false) * v2 -
            convert(Matrix{t}, EUR, center=true, model=model, impute=false) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, model=model, impute=false) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, scale=true, model=model, impute=false) * v2 -
            convert(Matrix{t}, EUR, scale=true, model=model, impute=false) * v2) /
            norm(convert(Matrix{t}, EUR, scale=true, model=model, impute=false) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=false) * v2 -
            convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false) * v2) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, model=model, impute=false)) * v1 -
            transpose(convert(Matrix{t}, EUR, model=model, impute=false)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, model=model, impute=false)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, center=true, model=model, impute=false)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=false)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=false)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, scale=true, model=model, impute=false)) * v1 -
            transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=false)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=false)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=false)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false)) * v1) < reltol
    end
end
end

@testset "SnpLinAlg-matrix multiplication zeroimpute" begin
    for t in [Float32, Float64] # Y = A*X
        X1 = randn(t, size(EUR, 1), 3)
        X2 = randn(t, size(EUR, 2), 3)
        for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
            @test all(isapprox(SnpLinAlg{t}(EUR, model=model, impute=false) * X2,
                convert(Matrix{t}, EUR, model=model, impute=false) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, center=true, model=model, impute=false) * X2,
                convert(Matrix{t}, EUR, center=true, model=model, impute=false) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, scale=true, model=model, impute=false) * X2,
                convert(Matrix{t}, EUR, scale=true, model=model, impute=false) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=false) * X2,
                convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false) * X2))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, model=model, impute=false)) * X1,
                transpose(convert(Matrix{t}, EUR, model=model, impute=false)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, center=true, model=model, impute=false)) * X1,
                transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=false)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, scale=true, model=model, impute=false)) * X1,
                transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=false)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=false)) * X1,
                transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=false)) * X1))
        end
    end
end

@testset "SnpLinAlg-vector multiplication (Miter > 0)" begin
    EUR11 = [EUR;EUR;EUR;EUR;EUR;EUR;EUR;EUR;EUR;EUR;EUR]
    EUR11la = SnpLinAlg{Float64}(EUR11, model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
    v = rand(size(EUR11la, 2))
    vtest = EUR11la * v
    vtrue = convert(Matrix{Float64}, EUR11, model=ADDITIVE_MODEL, impute=true, center=true, scale=true) * v
    @test norm(vtest - vtrue) < 5e-4
end

if get(ENV,"JULIA_SNPARRAYS_TEST_CUDA","") == "true"
    using CUDA, Adapt
    @testset "lin. alg. cuda" begin
    reltol = 5e-4
    for t in [Float32, Float64]
        v1 = randn(t, size(EUR, 1))
        v2 = randn(t, size(EUR, 2))
        v1_d = adapt(CuVector{t}, v1)
        v2_d = adapt(CuVector{t}, v2)
        for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
            @test norm(collect(CuSnpArray{t}(EUR, model=model) * v2_d) -
                convert(Matrix{t}, EUR, model=model) * v2) /
                norm(convert(Matrix{t}, EUR, model=model) * v2) < reltol
            @test norm(collect(CuSnpArray{t}(EUR, center=true, model=model) * v2_d) -
                convert(Matrix{t}, EUR, center=true, model=model) * v2) /
                norm(convert(Matrix{t}, EUR, center=true, model=model) * v2) < reltol
            @test norm(collect(CuSnpArray{t}(EUR, scale=true, model=model) * v2_d) -
                convert(Matrix{t}, EUR, scale=true, model=model) * v2) /
                norm(convert(Matrix{t}, EUR, scale=true, model=model) * v2) < reltol
            @test norm(collect(CuSnpArray{t}(EUR, center=true, scale=true, model=model) * v2_d) -
                convert(Matrix{t}, EUR, center=true, scale=true, model=model) * v2) /
                norm(convert(Matrix{t}, EUR, center=true, scale=true, model=model) * v2) < reltol
            @test norm(collect(transpose(CuSnpArray{t}(EUR, model=model)) * v1_d) -
                transpose(convert(Matrix{t}, EUR, model=model)) * v1) /
                norm(transpose(convert(Matrix{t}, EUR, model=model)) * v1) < reltol
            @test norm(collect(transpose(CuSnpArray{t}(EUR, center=true, model=model)) * v1_d) -
                transpose(convert(Matrix{t}, EUR, center=true, model=model)) * v1) /
                norm(transpose(convert(Matrix{t}, EUR, center=true, model=model)) * v1) < reltol
            @test norm(collect(transpose(CuSnpArray{t}(EUR, scale=true, model=model)) * v1_d) -
                transpose(convert(Matrix{t}, EUR, scale=true, model=model)) * v1) /
                norm(transpose(convert(Matrix{t}, EUR, scale=true, model=model)) * v1) < reltol
            @test norm(collect(transpose(CuSnpArray{t}(EUR, center=true, scale=true, model=model)) * v1_d) -
                transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model)) * v1) /
                norm(transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model)) * v1) < reltol
        end
    end
    end
end

@testset "SnpLinAlg-vector multiplication meanimpute" begin
reltol = 5e-4
for t in [Float32, Float64]
    v1 = randn(t, size(EUR, 1))
    v2 = randn(t, size(EUR, 2))
    for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
        @test norm(SnpLinAlg{t}(EUR, model=model, impute=true) * v2 -
            convert(Matrix{t}, EUR, model=model, impute=true) * v2) /
            norm(convert(Matrix{t}, EUR, model=model, impute=true) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, center=true, model=model, impute=true) * v2 -
            convert(Matrix{t}, EUR, center=true, model=model, impute=true) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, model=model, impute=true) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, scale=true, model=model, impute=true) * v2 -
            convert(Matrix{t}, EUR, scale=true, model=model, impute=true) * v2) /
            norm(convert(Matrix{t}, EUR, scale=true, model=model, impute=true) * v2) < reltol
        @test norm(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=true) * v2 -
            convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true) * v2) /
            norm(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true) * v2) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, model=model, impute=true)) * v1 -
            transpose(convert(Matrix{t}, EUR, model=model, impute=true)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, model=model, impute=true)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, center=true, model=model, impute=true)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=true)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=true)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, scale=true, model=model, impute=true)) * v1 -
            transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=true)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=true)) * v1) < reltol
        @test norm(transpose(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=true)) * v1 -
            transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true)) * v1) /
            norm(transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true)) * v1) < reltol
    end
end
end

@testset "SnpLinAlg-matrix multiplication meanimpute" begin
    for t in [Float32, Float64] # Y = A*X
        X1 = randn(t, size(EUR, 1), 3)
        X2 = randn(t, size(EUR, 2), 3)
        for model in [ADDITIVE_MODEL, DOMINANT_MODEL, RECESSIVE_MODEL]
            @test all(isapprox(SnpLinAlg{t}(EUR, model=model, impute=true) * X2,
                convert(Matrix{t}, EUR, model=model, impute=true) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, center=true, model=model, impute=true) * X2,
                convert(Matrix{t}, EUR, center=true, model=model, impute=true) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, scale=true, model=model, impute=true) * X2,
                convert(Matrix{t}, EUR, scale=true, model=model, impute=true) * X2))
            @test all(isapprox(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=false) * X2,
                convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true) * X2))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, model=model, impute=true)) * X1,
                transpose(convert(Matrix{t}, EUR, model=model, impute=true)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, center=true, model=model, impute=true)) * X1,
                transpose(convert(Matrix{t}, EUR, center=true, model=model, impute=true)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, scale=true, model=model, impute=true)) * X1,
                transpose(convert(Matrix{t}, EUR, scale=true, model=model, impute=true)) * X1))
            @test all(isapprox(transpose(SnpLinAlg{t}(EUR, center=true, scale=true, model=model, impute=true)) * X1,
                transpose(convert(Matrix{t}, EUR, center=true, scale=true, model=model, impute=true)) * X1))
        end
    end
end

@testset "subarrays" begin
@test all(@view(EUR[1:2:10, 1:2:10]) .==
[[0x03 0x03 0x02 0x02 0x03];
[0x03 0x03 0x03 0x03 0x03];
[0x03 0x03 0x03 0x03 0x03];
[0x02 0x03 0x02 0x00 0x02];
[0x03 0x03 0x02 0x02 0x03]
])
@test all(convert(Matrix{Float64}, @view(EUR[1:2:10, 1:2:10])) .≈
[[2.0 2.0 1.0 1.0 2.0];
[2.0 2.0 2.0 2.0 2.0];
[2.0 2.0 2.0 2.0 2.0];
[1.0 2.0 1.0 0.0 1.0];
[2.0 2.0 1.0 1.0 2.0]
])

EURsub = @view EUR[1:2:100, 1:2:100]
EURsubbm = SnpBitMatrix{Float64}(EURsub, model=ADDITIVE_MODEL, center=true, scale=true) # BitMatrix from the SubArray
EURsubfm = convert(Matrix{Float64}, EURsub, model=ADDITIVE_MODEL, center=true, scale=true) # Float64 Matrix from the SubArray

v1 = randn(size(EURsub, 1))
v2 = randn(size(EURsub, 2))

@test isapprox(EURsubbm * v2, EURsubfm * v2)
@test isapprox(EURsubbm' * v1, EURsubfm' * v1)
end

@testset "split-merge-readwrite" begin
EUR_data = SnpData(SnpArrays.datadir("EUR_subset"))
# a small subset of EUR_data
SnpArrays.filter(EUR_data, collect(1:10), collect(1:20000); des="EUR_subset.small")
EUR_small = SnpData("EUR_subset.small")

# filter
chr17 = SnpArrays.filter(SnpArrays.datadir("EUR_subset"); des="tmp.filter.chr17", f_snp = x -> String(x[:chromosome]) == "17")
@test chr17.snps == 11041
@test chr17.people == 379
male = SnpArrays.filter(EUR_data; des="tmp.filter.male", f_person = x -> String(x[:sex]) == "1")
@test male.snps == 54051
@test male.people == 178
chr17_male = SnpArrays.filter(EUR_data; des="tmp.filter.chr17.male", f_snp = x -> String(x[:chromosome]) == "17", f_person = x -> String(x[:sex]) == "1")
@test chr17_male.snps == 11041
@test chr17_male.people == 178
@test size(chr17_male.snparray) == (178, 11041)
# cleanup
for ft in ["bim", "fam", "bed"]
    ft == "bed" && Sys.iswindows() && continue
    rm("tmp.filter.chr17." * ft, force=true)
    rm("tmp.filter.male." * ft, force=true)
    rm("tmp.filter.chr17.male." * ft, force=true)
end

# split
splitted = SnpArrays.split_plink("EUR_subset.small"; prefix="tmp.split.chr.")
piece = splitted["17"]
@test piece.people == 10
@test piece.snps == 11041
@test size(piece.person_info) == (10, 6)
@test size(piece.snp_info) == (11041, 6)
@test size(piece.snparray) == (10, 11041)
@test size(piece.snparray.columncounts) == (4, 11041)
@test piece.snparray.m == 10

splitted_bysex = SnpArrays.split_plink(EUR_data, :sex; prefix="tmp.split.sex.")
@test splitted_bysex["1"].people == 178
@test splitted_bysex["2"].people == 201

# merge
#@time merged = SnpArrays.merge_plink("tmp.merged", splitted) # write_plink is included here
#@test EUR_data.people == merged.people
#@test EUR_data.snps == merged.snps
#@test EUR_data.person_info == merged.person_info
#@test EUR_data.snp_info == merged.snp_info # note: the ordering of merged data might be different on other dataset b/c sorted order of chromosomes

#output = SnpData("tmp.merged")
#@test EUR_data.people == output.people
#@test EUR_data.snps == output.snps
#@test EUR_data.person_info == output.person_info
#@test EUR_data.snp_info == output.snp_info

# cleanup
#rm("tmp.merged.bim", force=true)
#rm("tmp.merged.fam", force=true)
#Sys.iswindows() || rm("tmp.merged.bed", force=true)

# merge from splitted files
@time merged_from_splitted_files = merge_plink("tmp.split.chr"; des = "tmp2.merged")
@test EUR_small.people == merged_from_splitted_files.people
@test EUR_small.snps == merged_from_splitted_files.snps
@test EUR_small.person_info == merged_from_splitted_files.person_info
@test EUR_small.snp_info == merged_from_splitted_files.snp_info

# cleanup
rm("tmp2.merged.bim", force=true)
rm("tmp2.merged.fam", force=true)
Sys.iswindows() || rm("tmp2.merged.bed", force=true)

rm("EUR_subset.small.bim", force=true)
rm("EUR_subset.small.fam", force=true)
Sys.iswindows() || rm("EUR_subset.small.bed", force=true)

for k in keys(splitted)
    for ft in ["bim", "fam", "bed"]
        ft == "bed" && Sys.iswindows() && continue
        rm("tmp.split.chr.$(k)." * ft, force=true)
    end
end
for k in keys(splitted_bysex)
    for ft in ["bim", "fam", "bed"]
        ft == "bed" && Sys.iswindows() && continue
        rm("tmp.split.sex.$(k)." * ft, force=true)
    end
end
end

@testset "(de)compress" begin
for format in SnpArrays.ALLOWED_FORMAT
    # compress mouse Plink files
    compress_plink(SnpArrays.datadir("mouse"), format)
    @test isfile(SnpArrays.datadir("mouse.bed." * format))
    @test isfile(SnpArrays.datadir("mouse.fam." * format))
    @test isfile(SnpArrays.datadir("mouse.bim." * format))
    # read in compressed Plink files
    mouse_zip = SnpArray(SnpArrays.datadir("mouse.bed." * format))
    fill!(mouse.rowcounts, 0)
    fill!(mouse.columncounts, 0)
    @time begin
        SnpArrays._counts(mouse, 1)
        SnpArrays._counts(mouse, 2)
    end
    fill!(mouse_zip.rowcounts, 0)
    fill!(mouse_zip.columncounts, 0)
    @time begin
        SnpArrays._counts(mouse_zip, 1)
        SnpArrays._counts(mouse_zip, 2)
    end
    @test mouse.m == mouse_zip.m
    @test norm(mouse.rowcounts - mouse_zip.rowcounts) < 1e-8
    @test norm(mouse.columncounts - mouse_zip.columncounts) < 1e-8
    @test norm(mouse.data - mouse_zip.data) < 1e-8
    # Decompress Plink files
    decompress_plink(SnpArrays.datadir("mouse"), format, SnpArrays.datadir("mouse2"))
    @test stat(SnpArrays.datadir("mouse") * ".bed").size ==
    stat(SnpArrays.datadir("mouse2") * ".bed").size
    @test stat(SnpArrays.datadir("mouse") * ".fam").size ==
    stat(SnpArrays.datadir("mouse2") * ".fam").size
    @test stat(SnpArrays.datadir("mouse") * ".bim").size ==
    stat(SnpArrays.datadir("mouse2") * ".bim").size
    # clean up
    rm(SnpArrays.datadir("mouse.bed." * format), force=true)
    rm(SnpArrays.datadir("mouse.fam." * format), force=true)
    rm(SnpArrays.datadir("mouse.bim." * format), force=true)
    rm(SnpArrays.datadir("mouse2.bed"), force=true)
    rm(SnpArrays.datadir("mouse2.fam"), force=true)
    rm(SnpArrays.datadir("mouse2.bim"), force=true)
end
@test_throws ArgumentError SnpArray(SnpArrays.datadir("mouse.bed.zip"))
end

@testset "indexin" begin
    b = ['1', '2', '3', '4']
    for a in [['4', '3', '1', '2'], ['3', '2', '6'], ['b', 'c', 'd', 'c']]
        aind, bmask = SnpArrays.indexin_general(a, b)
        @test all(a[aind] .== b[bmask])
    end
end

@testset "concat" begin
s = SnpArrays.filter(SnpArrays.datadir("mouse"), 1:2, 1:3; des="mouse.tmp.filtered")
sd = SnpData("mouse.tmp.filtered")
sd_hcat = hcat(sd, sd, sd; des="mouse.hcat")
sd_vcat = vcat(sd, sd, sd; des="mouse.vcat")
sd_hvcat = hvcat((2, 2), sd, sd, sd, sd; des="mouse.hvcat")

ref = [[0x02 0x02 0x02];
[0x02 0x02 0x03]]

@test all(sd_hcat.snparray .== [ref ref ref])
@test all(sd_vcat.snparray .== [ref; ref; ref])
@test all(sd_hvcat.snparray .== [ref ref; ref ref])

@test sd_hcat.people == 2
@test sd_hcat.snps == 9

@test sd_vcat.people == 6
@test sd_vcat.snps == 3

@test sd_hvcat.people == 4
@test sd_hvcat.snps == 6

for ft in ["bim", "fam", "bed"]
    ft == "bed" && Sys.iswindows() && continue
    rm("mouse.tmp.filtered." * ft, force=true)
    rm("mouse.hcat." * ft, force=true)
    rm("mouse.vcat." * ft, force=true)
    rm("mouse.hvcat." * ft, force=true)
end
end

@testset "reorder" begin
mouse_prefix = SnpArrays.datadir("mouse")
run(`cp $(mouse_prefix * ".bed") mouse_testreorder.bed`)
run(`cp $(mouse_prefix * ".bim") mouse_testreorder.bim`)
run(`cp $(mouse_prefix * ".fam") mouse_testreorder.fam`)

mouse_data = SnpData(mouse_prefix)
mouse_toreorder = SnpData("mouse_testreorder", "r+")
m, n = size(mouse_toreorder.snparray)
using Random
ind = randperm(m);
SnpArrays.reorder!(mouse_toreorder, ind)

# reread the file
mouse_toreorder_read = SnpData("mouse_testreorder"; famnm="mouse_testreorder.reordered.fam")
@test all(mouse_data.snparray[ind, :] .== mouse_toreorder_read.snparray[:, :])
@test all(mouse_data.person_info[ind, 1] .== mouse_toreorder_read.person_info[:, 1])

# cleanup
rm("mouse_testreorder.bed", force=true)
rm("mouse_testreorder.bim", force=true)
rm("mouse_testreorder.fam", force=true)
rm("mouse_testreorder.reordered.fam", force=true)
end

@testset "vcf2plink" begin
# Download an example VCF file
isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
    joinpath(pwd(), "test.08Jun17.d8b.vcf.gz"));
s = vcf2plink("test.08Jun17.d8b.vcf.gz", "test.08Jun17.d8b")
@test size(s, 1) == 191
@test size(s, 2) == 1354
@test s[1, 1] == 0x00
@test s[1, 5] == 0x02
@test s[38, 5] == 0x03
rm("test.08Jun17.d8b.vcf.gz", force=true)
rm("test.08Jun17.d8b.bed", force=true)
rm("test.08Jun17.d8b.bim", force=true)
rm("test.08Jun17.d8b.fam", force=true)
end

@testset "kinship_pruning" begin
g = grm(mouse)
@test count(kinship_pruning(g; method=:gcta)) == 68
@test count(kinship_pruning(g; method=:top_down)) == 125
@test count(kinship_pruning(g; method=:bottom_up)) == 132
@test count(kinship_pruning(g; method=:plink)) == 126
end

@testset "stackedsnparray" begin
mouse2 = StackedSnpArray([mouse, mouse])
@test size(mouse2) == (1940, 20300)
@test size(mouse2, 3) == 1
@test_throws ErrorException size(mouse2, 0)
@test length(mouse2) == 1940 * 20300
@test eltype(mouse2) == UInt8
@test mouse2[777, 16384] == mouse[777, 16384 - 10150]
end

@testset "SNP simulation" begin
    function observed_mafs!(
        obs::Vector{T},
        counts::Matrix{Int}
        ) where {T}
        @inbounds for j in axes(counts, 2)
            obs[j] = (counts[3, j] + 2 * counts[4, j]) / (2 * (counts[1, j] + counts[3, j] + counts[4, j]))
        end
    end
    rng = MersenneTwister(1234)
    reltol = 1e-2
    M = 1_024 * 1000
    N = 1_024
    minfreq = 0.01
    MAFs = minfreq .+ (0.5 - minfreq) * rand(rng, N)
    G = SnpArrays.simulate(rng, M, N, MAFs)
    # Simulate missing values
    missing_rates = 0.10 * rand(rng, N)
    SnpArrays.simulate_missing!(rng, G, missing_rates)
    # Get counts
    cc = counts(G; dims=1)
    observed_prob = similar(MAFs);
    observed_mafs!(observed_prob, cc)
    @test norm(observed_prob - MAFs) / norm(MAFs) < reltol
    observed_missing_rates = cc[2, :] ./ M
    @test norm(observed_missing_rates - missing_rates) / norm(missing_rates) < reltol
    # Simulation under linkage disequilibrium
    ρ = 0.10
    G = SnpArrays.simulate!(rng, G, MAFs, ρ)
    cc = counts(G; dims=1)
    observed_mafs!(observed_prob, cc)
    @test norm(observed_prob - MAFs) / norm(MAFs) < reltol
    # check for first order correlation
    σ = zeros(N)
    for i in 1:N
        σ[i] = sqrt(2 * observed_prob[i] * (1 - observed_prob[i]))
    end
    r = zeros(N - 1)
    g1, g2 = zeros(M), zeros(M)
    for j in eachindex(r)
        copyto!(g1, view(G, :, j))
        g1 .-= 2 * observed_prob[j]
        copyto!(g2, view(G, :, j + 1))
        g2 .-= 2 * observed_prob[j + 1]
        r[j] = dot(g1, g2) / ((M-1) * σ[j] * σ[j + 1])
    end
    @test abs(mean(r) - ρ) / ρ < reltol
end
