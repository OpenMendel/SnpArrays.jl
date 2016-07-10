module HaplotypeArraysTest

using SnpArrays
if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

info("Test HaplotypeArray implementation")

srand(123)
n, p = 100, 1000

@testset "Constructors" begin
  a1 = rand(0:1, n, p)
  a2 = rand(0:1, n, p)
  h = HaplotypeArray(a1 .> 0, a2 .> 0)
  @test typeof(h) <: AbstractMatrix
  @test size(h) == size(a1)
  @test eltype(h) == Tuple{Bool, Bool}
  @test endof(h) == n * p
  @test ndims(h) == 2
  # construct SnpArray from HaplotypeArray
  snp = SnpArray(h)
  @test any(isnan(snp)) == false
  # construct HaplotypeArray from SnpArray
  h2 = HaplotypeArray(snp)
  @test h2 == snp # bitwise same
end

@testset "{g,s}etindex" begin
  snp = HaplotypeArray(n, p)
  i, j = rand(1:n), rand(1:p)
  @test snp[i, j] == snp[sub2ind((n, p), i, j)]
  a1, a2 = rand(Bool), rand(Bool)
  snp[i, j] = (a1, a2)
  @test snp[i, j] == (a1, a2)
end

@testset "convert" begin
  #@code_llvm convert(Float64, (false, false), true, :additive)
  # A1 is the minor allele
  @test convert_haplotype(Float64, (false, false), true, :additive) == 2.0
  @test convert_haplotype(Float64, (false, true), true, :additive) == 1.0
  @test convert_haplotype(Float64, (true, false), true, :additive) == 1.0
  @test convert_haplotype(Float64, (true, true), true, :additive) == 0.0
  @test convert_haplotype(Float64, (false, false), true, :dominant) == 1.0
  @test convert_haplotype(Float64, (false, true), true, :dominant) == 1.0
  @test convert_haplotype(Float64, (true, false), true, :dominant) == 1.0
  @test convert_haplotype(Float64, (true, true), true, :dominant) == 0.0
  @test convert_haplotype(Float64, (false, false), true, :recessive) == 1.0
  @test convert_haplotype(Float64, (false, true), true, :recessive) == 0.0
  @test convert_haplotype(Float64, (true, false), true, :recessive) == 0.0
  @test convert_haplotype(Float64, (true, true), true, :recessive) == 0.0
  # A2 is the minor allele
  @test convert_haplotype(Float64, (false, false), false, :additive) == 0.0
  @test convert_haplotype(Float64, (false, true), false, :additive) == 1.0
  @test convert_haplotype(Float64, (true, false), true, :additive) == 1.0
  @test convert_haplotype(Float64, (true, true), false, :additive) == 2.0
  @test convert_haplotype(Float64, (false, false), false, :dominant) == 0.0
  @test convert_haplotype(Float64, (false, true), false, :dominant) == 1.0
  @test convert_haplotype(Float64, (true, false), false, :dominant) == 1.0
  @test convert_haplotype(Float64, (true, true), false, :dominant) == 1.0
  @test convert_haplotype(Float64, (false, false), false, :recessive) == 0.0
  @test convert_haplotype(Float64, (false, true), false, :recessive) == 0.0
  @test convert_haplotype(Float64, (true, false), false, :recessive) == 0.0
  @test convert_haplotype(Float64, (true, true), false, :recessive) == 1.0
  # convert to real matrix
  h = HaplotypeArray(n, p)
  h_float = convert(Matrix{Float64}, h)
  @test typeof(h_float) == Matrix{Float64}
  @test size(h) == size(h_float)
  # copy to real matrix
  storage = zeros(n, p)
  #@code_warntype copy!(storage, snp)
  copy!(storage, h)
  @test all(h_float .== storage)
  # convert to sparse matrix
  h_sparse = SparseMatrixCSC{Float64, Int}(h)
  @test typeof(h_sparse) == SparseMatrixCSC{Float64, Int}
  @test h_sparse == h_float
end

@testset "summarize" begin
  # without missing genotypes
  maf, minor_allele = 0.5rand(p), bitrand(p)
  snp = randgeno(n, p, maf, minor_allele)
  h = HaplotypeArray(snp)
  #@code_warntype summarize(snp)
  @inferred summarize(h)
  maf, minor_allele = summarize(h)
  @test all(0.0 .≤ maf .≤ 1.0)
  @test eltype(minor_allele) == Bool
  # summarize a HaplotypeVector
  maf, minor_allele = summarize(h[:, 1])
  @test 0.0 ≤ maf ≤ 1.0
  @test typeof(minor_allele) == Bool
  # summarize a view of HaplotypeVector
  h1 = h[:, 1]
  summarize(sub(h1, 1:10))
  # corner case: m = 0 (SnpArray is empty)
  h = HaplotypeArray(0, 5)
  maf, = summarize(h)
  @test isnan(maf[1])
  maf, = summarize(h[:, 1])
  @test isnan(maf[1])
  # edge case: p = 0 (SnpArray is empty)
  h = HaplotypeArray(5, 0)
  maf, = summarize(h)
  @test isempty(maf)
end


end # HaplotypeArraysTest
