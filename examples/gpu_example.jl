using OpenCL
using SnpArrays
using LinearAlgebra, BenchmarkTools

# setting up
const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"));
using Random
Random.seed!(95376)
X = EUR
m,n  = size(EUR)

X_bm = SnpBitMatrix{Float64}(X; model=ADDITIVE_MODEL, center=true, scale=true);
y = randn(m)
z = Vector{Float64}(undef, n)
z_bm = Vector{Float64}(undef, n)
y32 = convert(Array{Float32}, y)
z32 = convert(Array{Float32}, z)

v = SnpCLVariables(z, X, y)
v32 = SnpCLVariables(z32, X, y32)

# benchmarks
print("Benchmark time for CPU version (via SnpBitMatrix): ")
mul!(z_bm, transpose(X_bm), y)
@btime mul!(z_bm, transpose(X_bm), y)

print("Benchmark time for GPU version (64-bit): ")
mul!(z, transpose(X), y; v=v)
@btime mul!(z, transpose(X), y; v=v)

print("Benchmark time for GPU version (32-bit): ")
mul!(z32, transpose(X), y32; v=v32)
@btime mul!(z32, transpose(X), y32; v=v32)

# check correctness
@assert isapprox(z, z_bm; rtol=1e-12)
@assert isapprox(z32, z_bm; rtol=1e-5)
@info "SUCCESS!"
