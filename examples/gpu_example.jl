using OpenCL
using SnpArrays
using LinearAlgebra, BenchmarkTools, Printf

# platform an device information
platform = first(cl.platforms())
@printf("Platform name:    %s\n",  platform[:name])
@printf("Platform profile: %s\n",  platform[:profile])
@printf("Platform vendor:  %s\n",  platform[:vendor])
@printf("Platform version: %s\n",  platform[:version])

println()
device = last(cl.devices(:gpu))
@printf("Device name: %s\n", device[:name])
@printf("Device type: %s\n", device[:device_type])
@printf("Device mem: %i MB\n",           device[:global_mem_size] / 1024^2)
@printf("Device max mem alloc: %i MB\n", device[:max_mem_alloc_size] / 1024^2)
@printf("Device max clock freq: %i MHZ\n",  device[:max_clock_frequency])
@printf("Device max compute units: %i\n",   device[:max_compute_units])
@printf("Device max work group size: %i\n", device[:max_work_group_size])
@printf("Device max work item size: %s\n",  device[:max_work_item_size]) 


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
