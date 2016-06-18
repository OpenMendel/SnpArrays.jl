module SnpArraysTest
using SnpArrays
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

srand(123)
n, p = 100, 1000

@testset "Constructors" begin
  x = rand(0:2, n, p)
  snp = SnpArray(x)
  @test typeof(snp) <: AbstractMatrix
  @test size(snp) == size(x)
  @test eltype(snp) == Tuple{Bool, Bool}
end

@testset "{g,s}etindex" begin
  snp = SnpArray(n, p)
  i, j = rand(1:n), rand(1:p)
  snp[i, j] = 2
  @test snp[i, j] == (true, true)
  snp[i, j] = 1
  @test snp[i, j] == (false, true)
  snp[i, j] = 0
  @test snp[i, j] == (false, false)
  snp[i, j] = NaN
  @test snp[i, j] == (true, false)
end

@testset "convert" begin
  snp = SnpArray(n, p)
  snp_float = Matrix{Float64}(snp)
  @test typeof(snp_float) == Matrix{Float64}
  @test size(snp) == size(snp_float)
  storage = zeros(n, p)
  copy!(storage, snp)
  @test snp_float == storage
  snp_sparse = SparseMatrixCSC{Float64, Int}(snp)
  @test typeof(snp_sparse) == SparseMatrixCSC{Float64, Int}
  @test snp_sparse == snp_float
end

@testset "summarize" begin
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  _, _, missings_by_snp, missings_by_person = summarize(snp)
  @test sum(missings_by_snp) == 0
  @test sum(missings_by_person) == 0
  @test sum(isnan(snp)) == 0
end

@testset "grm" begin
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  # empirical kinship by GRM
  Φ_grm = grm(snp; method = :GRM)
  @test issym(Φ_grm)
  @test all(eigvals(Φ_grm) .≥ -1.0e-6)
  # empirical kinship by MoM
  Φ_mom = grm(snp; method = :MoM)
  @test issym(Φ_mom)
  # MoM is not guaranteed to be psd:
  #@test all(eigvals(Φ_mom) .≥ -1.0e-6)
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
  @test vecnorm(pcvariance_f32 - pcvariance_f64) / vecnorm(pcvariance_f64) ≤ 1.0e-4
end

end # module
