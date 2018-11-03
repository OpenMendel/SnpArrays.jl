using SnpArrays, SparseArrays, LinearAlgebra, Test

const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"))
const mouse = SnpArray(SnpArrays.datadir("mouse.bed"))

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
    rm("tmp.bed")
    tmpbf = SnpArray("tmp.bed", SnpArray(undef, 5, 3))
    fill!(tmpbf, 0x01)
    @test all(tmpbf .== 0x01)
    rm("tmp.bed")
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
    @test all(eigvals(Φgrm) .≥ 0)
    Φmom = grm(EUR, method=:MoM)
    @test size(Φmom) == (size(EUR, 1), size(EUR, 1))
    @test issymmetric(Φmom)
    Φrbs = grm(EUR, method=:Robust)
    @test size(Φrbs) == (size(EUR, 1), size(EUR, 1))
    @test issymmetric(Φrbs)
    @test all(eigvals(Φrbs) .≥ 0)
end

@testset "filter" begin
    rowmask, colmask =  SnpArrays.filter(mouse, 0.99, 0.99)
    SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask; des="tmp")
    tmpbf = SnpArray("tmp.bed")
    @test size(tmpbf) == (1907, 9997)
    @test all(missingrate(tmpbf, 1) .≤ 0.01)
    @test all(missingrate(tmpbf, 2) .≤ 0.01)
    rm("tmp.bed")
    rm("tmp.bim")
    rm("tmp.fam")
end
