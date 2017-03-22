module SnpArraysTest

using Compat
import Compat:view, is_unix, issymmetric

using SnpArrays
if VERSION >= v"0.5.0"
    using Base.Test
else
    using BaseTestNext
end

info("Test SnpArray implementation")

srand(123)
n, p = 100, 1000

# @testset "Constructors" begin
#   x = rand(0:2, n, p)
#   snp = SnpArray(x)
#   @test typeof(snp) <: AbstractMatrix
#   @test size(snp) == size(x)
#   @test eltype(snp) == Tuple{Bool, Bool}
#   @test endof(snp) == n * p
#   # constructor from Plink files
#   hapmap1 = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3")
#   @test size(hapmap1) == (324, 13928)
#   hapmap2 = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3"; people = 324, snps = 13928)
#   @test all(hapmap1 == hapmap2)
# end

@testset "{g,s}etindex" begin
  snp = SnpArray(n, p)
  i, j = rand(1:n), rand(1:p)
  @test snp[i, j] == snp[sub2ind((n, p), i, j)]
  snp[i, j] = 2
  @test snp[i, j] == (true, true)
  snp[i, j] = 2.0
  @test snp[i, j] == (true, true)
  snp[i, j] = 1
  @test snp[i, j] == (false, true)
  snp[i, j] = 1.0
  @test snp[i, j] == (false, true)
  snp[i, j] = 0
  @test snp[i, j] == (false, false)
  snp[i, j] = 0.0
  @test snp[i, j] == (false, false)
  snp[i, j] = NaN
  @test snp[i, j] == (true, false)
  snp[sub2ind((n, p), i, j)] = 2.0
  @test snp[sub2ind((n, p), i, j)] == (true, true)
  snp[sub2ind((n, p), i, j)] = NaN
  @test snp[sub2ind((n, p), i, j)] == (true, false)
end

@testset "similar" begin
  snp = SnpArray(n, p)
  A = similar(snp, Float64, size(snp))
  @test eltype(A) == Float64
  @test size(A) == size(snp)
  A = similar(snp)
  @test typeof(A) == SnpArray{2}
  @test size(A) == size(snp)
  A = similar(snp, (3, 3))
  @test typeof(A) == SnpArray{2}
  @test size(A) == (3, 3)
  A = similar(snp, Tuple{Bool, Bool})
  @test typeof(A) == SnpArray{2}
  @test size(A) == size(snp)
end

@testset "isnan" begin
  snp = SnpArray(n, p)
  #@code_warntype isnan(snp)
  @inferred isnan(snp)
  @test isnan((true, false)) == true
  @test isnan((true, true)) == false
  @test isnan((false, true)) == false
  @test isnan((false, false)) == false
end

@testset "convert" begin
  #@code_llvm convert(Float64, (false, false), true, :additive)
  # A1 is the minor allele
  @test convert(Float64, (false, false), true, :additive) == 2.0
  @test convert(Float64, (false, true), true, :additive) == 1.0
  @test convert(Float64, (true, true), true, :additive) == 0.0
  @test isnan(convert(Float64, (true, false), true, :additive))
  @test convert(Float64, (false, false), true, :dominant) == 1.0
  @test convert(Float64, (false, true), true, :dominant) == 1.0
  @test convert(Float64, (true, true), true, :dominant) == 0.0
  @test isnan(convert(Float64, (true, false), true, :dominant))
  @test convert(Float64, (false, false), true, :recessive) == 1.0
  @test convert(Float64, (false, true), true, :recessive) == 0.0
  @test convert(Float64, (true, true), true, :recessive) == 0.0
  @test isnan(convert(Float64, (true, false), true, :recessive))
  # A2 is the minor allele
  @test convert(Float64, (false, false), false, :additive) == 0.0
  @test convert(Float64, (false, true), false, :additive) == 1.0
  @test convert(Float64, (true, true), false, :additive) == 2.0
  @test isnan(convert(Float64, (true, false), true, :additive))
  @test convert(Float64, (false, false), false, :dominant) == 0.0
  @test convert(Float64, (false, true), false, :dominant) == 1.0
  @test convert(Float64, (true, true), false, :dominant) == 1.0
  @test isnan(convert(Float64, (true, false), false, :dominant))
  @test convert(Float64, (false, false), false, :recessive) == 0.0
  @test convert(Float64, (false, true), false, :recessive) == 0.0
  @test convert(Float64, (true, true), false, :recessive) == 1.0
  @test isnan(convert(Float64, (true, false), false, :recessive))
  # convert to real matrix
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  snp_float = Matrix{Float64}(snp)
  @test typeof(snp_float) == Matrix{Float64}
  @test size(snp) == size(snp_float)
  @test estimatesize(n, p, Matrix{Float64}, mean(maf)) ≈ Base.summarysize(snp_float)
  # copy to real matrix
  storage = zeros(n, p)
  #@code_warntype copy!(storage, snp)
  copy!(storage, snp)
  @test all(snp_float .== storage)
  # convert to sparse matrix
  snp_sparse = SparseMatrixCSC{Float64, Int}(snp)
  @test typeof(snp_sparse) == SparseMatrixCSC{Float64, Int}
  @test snp_sparse == snp_float
  estimatesize(n, p, typeof(snp_sparse)), Base.summarysize(snp_sparse)
  # convert in presence of missing genotypes
  hapmap = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3")
  @test any(isnan(hapmap))
  @test any(isnan(convert(Matrix{Float64}, hapmap))) == true
  @test any(isnan(convert(Matrix{Float64}, hapmap; impute = true))) == false
  @test any(isnan(convert(SparseMatrixCSC{Float64, Int}, hapmap; impute = true))) == false
end

@testset "randgeno" begin
  @test typeof(randgeno(0.25, true)) == Tuple{Bool, Bool}
  @test typeof(randgeno(5, 0.25, true)) == SnpArray{1}
  @test typeof(randgeno((5, 5), 0.25ones(5), trues(5))) == SnpArray{2}
end

@testset "summarize" begin
  # without missing genotypes
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  #@code_warntype summarize(snp)
  @inferred summarize(snp)
  _, _, missings_by_snp, missings_by_person = summarize(snp)
  @test sum(missings_by_snp) == 0
  @test sum(missings_by_person) == 0
  @test sum(isnan(snp)) == 0
  # with missing genotypes
  snp[:, 1] = (true, false) # set first column to be all missing genotypes
  _, _, missings_by_snp, missings_by_person = summarize(snp)
  @test sum(missings_by_snp) == n
  @test sum(missings_by_person) == n
  # summarize a SnpVector
  _, _, missings = summarize(snp[:, 1])
  @test missings == n
  _, _, missings = summarize(view(snp, :, 1))
  @test missings == n
  # summarize a view of SnpVector
    snp1 = snp[:, 1]
  summarize(view(snp1, 1:10))
  # corner case: m = 0 (SnpArray is empty)
  snp = SnpArray(0, 5)
  maf, = summarize(snp)
  @test isnan(maf[1])
  maf, = summarize(snp[:, 1])
  @test isnan(maf[1])
  # edge case: p = 0 (SnpArray is empty)
  snp = SnpArray(5, 0)
  maf, = summarize(snp)
  @test isempty(maf)
end

@testset "filter SnpArray" begin
  hapmap = SnpArray(Pkg.dir("SnpArrays") * "/docs/hapmap3")
  #@code_warntype filter(hapmap, 0.98, 0.98, 5)
  @inferred filter(hapmap, 0.98, 0.98, 5)
  snp_idx, person_idx = filter(hapmap, 0.98, 0.98)
  _, _, missings_by_snp, missings_by_person =
    summarize(view(hapmap, person_idx, snp_idx))
  @test 1.0 - maximum(missings_by_snp) / length(missings_by_person) ≥ 0.98
  @test 1.0 - maximum(missings_by_person) / length(missings_by_snp) ≥ 0.98

  snpdata = SnpData(Pkg.dir("SnpArrays") * "/docs/hapmap3")
  filtered_data = filter(snpdata, snp_idx, person_idx)
  n, p = sum(person_idx), sum(snp_idx)
  @test filtered_data.people == n
  @test filtered_data.snps == p
  @test length(filtered_data.personid) == n
  @test length(filtered_data.snpid) == p
  @test length(filtered_data.chromosome) == p
  @test length(filtered_data.genetic_distance) == p
  @test length(filtered_data.basepairs) == p
  @test length(filtered_data.allele1) == p
  @test length(filtered_data.allele2) == p
  @test length(filtered_data.maf) == p
  @test length(filtered_data.minor_allele) == p
  @test size(filtered_data.snpmatrix) == (n, p)
  @test length(filtered_data.missings_per_person) == n
  @test length(filtered_data.missings_per_snp) == p
 end

@testset "filter SnpData" begin
  snpdata = SnpData(Pkg.dir("SnpArrays") * "/docs/hapmap3")
  filtered_data = filter(snpdata, 0.98, 0.98, 5)
  _, _, missings_by_snp, missings_by_person = summarize(filtered_data.snpmatrix)
  @test 1.0 - maximum(missings_by_snp) / length(missings_by_person) ≥ 0.98
  @test 1.0 - maximum(missings_by_person) / length(missings_by_snp) ≥ 0.98
end

@testset "grm" begin
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  # empirical kinship by GRM
  #@code_warntype _grm(snp, 2.0^30)
  @inferred _grm(snp, 2.0^30)
  Φ_grm = grm(snp; method = :GRM)
  @test @compat issymmetric(Φ_grm)
  @test all(eigvals(Φ_grm) .≥ -1.0e-6)
  # GRM: restrict memory usage to 5KB
  Φ_grm_memlim = grm(snp; method = :GRM, memory_limit = 5 * 2^10)
  @test_approx_eq_eps vecnorm(Φ_grm - Φ_grm_memlim) 0.0 1.0e-8
  # empirical kinship by MoM
  #@code_warntype _mom(snp, 2.0^30)
  @inferred _mom(snp, 2.0^30)
  Φ_mom = grm(snp; method = :MoM)
  @test @compat issymmetric(Φ_mom)
  # MoM is not guaranteed to be psd:
  #@test all(eigvals(Φ_mom) .≥ -1.0e-6)
  # MoM: restrict memory usage to 1KB
  Φ_mom_memlim = grm(snp; method = :MoM, memory_limit = 5 * 2^10)
  @test_approx_eq_eps vecnorm(Φ_mom - Φ_mom_memlim) 0.0 1.0e-8
end

@testset "pca" begin
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  # PCA
  _, _, pcvariance_f64 = pca(snp, 3)
  _, _, pcvariance_f32 = pca(snp, 3, Matrix{Float32})
  @test vecnorm(pcvariance_f32 - pcvariance_f64) / vecnorm(pcvariance_f64) ≤ 1.0e-4
  _, _, pcvariance_f16 = pca(snp, 3, Matrix{Float32})
  @test vecnorm(pcvariance_f16 - pcvariance_f64) / vecnorm(pcvariance_f64) ≤ 1.0e-4
  _, _, pcvariance_f32sp = pca_sp(snp, 3, SparseMatrixCSC{Float32, UInt32})
  @test vecnorm(pcvariance_f32sp - pcvariance_f64) / vecnorm(pcvariance_f64) ≤ 1.0e-4
end

@testset "SnpData" begin
  hapmapdata = SnpData(Pkg.dir("SnpArrays") * "/docs/hapmap3")
  @test size(hapmapdata.snpmatrix) == (324, 13928)
  @test hapmapdata.people == 324
  @test hapmapdata.snps == 13928
  @test length(hapmapdata.personid) == 324
  @test length(hapmapdata.snpid) == 13928
  @test length(hapmapdata.chromosome) == 13928
  @test length(hapmapdata.genetic_distance) == 13928
  @test length(hapmapdata.basepairs) == 13928
  @test length(hapmapdata.allele1) == 13928
  @test length(hapmapdata.allele2) == 13928
  @test length(hapmapdata.maf) == 13928
  @test length(hapmapdata.minor_allele) == 13928
  @test length(hapmapdata.missings_per_person) == 324
  @test length(hapmapdata.missings_per_snp) == 13928
end

# somehow the rm function causes error on Windows
@static if is_unix()
  @testset "Write Plink" begin
    hapmapdata = SnpData(Pkg.dir("SnpArrays") * "/docs/hapmap3")
    testfile = Pkg.dir("SnpArrays") * "/test/testhapmap"
    writeplink(testfile, hapmapdata)
    hapmaptest = SnpArray(testfile; people = 324, snps = 13928)
    @test all(hapmapdata.snpmatrix .== hapmaptest)
    rm(testfile * ".bed")
    rm(testfile * ".bim")
  end
end

end # SnpArraysTest module
