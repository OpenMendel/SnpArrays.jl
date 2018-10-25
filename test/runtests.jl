using BEDFiles, SparseArrays, Test

const EUR = BEDFile(BEDFiles.datadir("EUR_subset.bed"))
const mouse = BEDFile(BEDFiles.datadir("mouse.bed"))

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