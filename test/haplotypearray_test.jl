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
  snp = HaplotypeArray(a1 .> 0, a2 .> 0)
  @test typeof(snp) <: AbstractMatrix
  @test size(snp) == size(a1)
  @test eltype(snp) == Tuple{Bool, Bool}
end

end # HaplotypeArraysTest
