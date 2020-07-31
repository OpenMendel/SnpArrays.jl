# Linear algebra of SnpArray

SnpArrays.jl supports three modes of matrix-vector multiplications.

1. Direct operations on a plink-formatted `SnpArray`: `SnpLinAlg`
2. Operations on transformed `BitMatrix`es: `SnpBitMatrix`
3. Direct operations on a plink-formatted data on an Nvidia GPU: `CuSnpArray`.

`SnpLinAlg` and `SnpBitMatrix` use Chris Elrod's [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl) internally. It is much faster on machines with AVX support. `CuSnpArray` uses [CUDA.jl](https://juliagpu.gitlab.io/CUDA.jl/) internally.
On this page, we compare these three.


```julia
versioninfo()
```

    Julia Version 1.4.1
    Commit 381693d3df* (2020-04-14 17:20 UTC)
    Platform Info:
      OS: Linux (x86_64-pc-linux-gnu)
      CPU: Intel(R) Xeon(R) Silver 4114 CPU @ 2.20GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-8.0.1 (ORCJIT, skylake)



```julia
using SnpArrays
```


```julia
const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"));
```

Let's try with EUR data repeated 100 and 101 times: 37900 by 54051 and 38279 by 54051, respectively.


```julia
EUR_10 = [EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR]
EUR_100 = [EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10];
EUR_101 = [EUR_100; EUR];
```

We create instnaces of SnpLinAlg, SnpBitmatrix and CuSnpArray:


```julia
EUR_100_bm = SnpBitMatrix{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_100_sla = SnpLinAlg{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_100_sla_ = SnpLinAlg{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false, impute=false)
EUR_100_mat = convert(Matrix{Float64}, EUR_100, model=ADDITIVE_MODEL, center=false, scale=false);

EUR_101_bm = SnpBitMatrix{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_101_sla = SnpLinAlg{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_101_sla_ = SnpLinAlg{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false, impute=false)
EUR_101_mat = convert(Matrix{Float64}, EUR_101, model=ADDITIVE_MODEL, center=false, scale=false);
```


```julia
ENV["JULIA_CUDA_USE_BINARYBUILDER"] = "false"
using CUDA
EUR_100_cu = CuSnpArray{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_100_cu_ = CuSnpArray{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false, impute=false);
```

    ┌ Warning: `haskey(::TargetIterator, name::String)` is deprecated, use `Target(; name = name) !== nothing` instead.
    │   caller = llvm_compat(::VersionNumber) at compatibility.jl:176
    └ @ CUDA /home/kose/.julia/packages/CUDA/5t6R9/deps/compatibility.jl:176


## $y = Ax$


```julia
using LinearAlgebra
using BenchmarkTools
```


```julia
v1 = randn(size(EUR_100, 1))
v1_ = randn(size(EUR_100, 1))
v2 = randn(size(EUR_100, 2));
```

With 8-threaded OpenBLAS (included in standard binary installation of Julia): 


```julia
BLAS.set_num_threads(8)
@benchmark LinearAlgebra.mul!($v1, $EUR_100_mat, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  0 bytes
      allocs estimate:  0
      --------------
      minimum time:     383.853 ms (0.00% GC)
      median time:      417.986 ms (0.00% GC)
      mean time:        431.205 ms (0.00% GC)
      maximum time:     530.847 ms (0.00% GC)
      --------------
      samples:          12
      evals/sample:     1



With single-threaded OpenBLAS: 


```julia
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($v1, $EUR_100_mat, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  0 bytes
      allocs estimate:  0
      --------------
      minimum time:     1.849 s (0.00% GC)
      median time:      2.064 s (0.00% GC)
      mean time:        2.047 s (0.00% GC)
      maximum time:     2.229 s (0.00% GC)
      --------------
      samples:          3
      evals/sample:     1



Direct linear algebra on a SnpArray, with mean imputation: 


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_100_sla, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  38.33 KiB
      allocs estimate:  1616
      --------------
      minimum time:     1.241 s (0.00% GC)
      median time:      1.355 s (0.00% GC)
      mean time:        1.352 s (0.00% GC)
      maximum time:     1.457 s (0.00% GC)
      --------------
      samples:          4
      evals/sample:     1



With zero imputation:


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_100_sla_, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  38.33 KiB
      allocs estimate:  1616
      --------------
      minimum time:     1.017 s (0.00% GC)
      median time:      1.037 s (0.00% GC)
      mean time:        1.034 s (0.00% GC)
      maximum time:     1.045 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1



Indeed, we are paying some price for mean imputation.

The below is the benchmark for SnpBitMatrix (always zero-imputed):


```julia
@benchmark (LinearAlgebra.mul!($v1, $EUR_100_bm, $v2))
```




    BenchmarkTools.Trial: 
      memory estimate:  0 bytes
      allocs estimate:  0
      --------------
      minimum time:     1.045 s (0.00% GC)
      median time:      1.051 s (0.00% GC)
      mean time:        1.051 s (0.00% GC)
      maximum time:     1.058 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1



At first glance, the result from SnpBitMatrix might look better than SnpLinAlg. However, SnpLinAlg is more stable in performance when the number of samples is not multiple of 4 or 8.


```julia
v1 = randn(size(EUR_101, 1))
v2 = randn(size(EUR_101, 2));
```


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_sla, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  44.13 KiB
      allocs estimate:  1722
      --------------
      minimum time:     1.259 s (0.00% GC)
      median time:      1.268 s (0.00% GC)
      mean time:        1.268 s (0.00% GC)
      maximum time:     1.277 s (0.00% GC)
      --------------
      samples:          4
      evals/sample:     1




```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_sla_, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  44.13 KiB
      allocs estimate:  1722
      --------------
      minimum time:     1.038 s (0.00% GC)
      median time:      1.043 s (0.00% GC)
      mean time:        1.043 s (0.00% GC)
      maximum time:     1.051 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1




```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_bm, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  0 bytes
      allocs estimate:  0
      --------------
      minimum time:     1.544 s (0.00% GC)
      median time:      1.561 s (0.00% GC)
      mean time:        1.561 s (0.00% GC)
      maximum time:     1.577 s (0.00% GC)
      --------------
      samples:          4
      evals/sample:     1



Now let's try CUDA. The device is Nvidia Titan V.


```julia
using Adapt
```

Moving data to GPU: 


```julia
v1 = randn(size(EUR_100, 1))
v1_ = randn(size(EUR_100, 1))
v2 = randn(size(EUR_100, 2));
v1_d = adapt(CuArray{Float64}, v1)
v1_d_ = similar(v1_d)
v2_d = adapt(CuArray{Float64}, v2);
```


```julia
using BenchmarkTools
@benchmark CUDA.@sync LinearAlgebra.mul!($v1_d, $EUR_100_cu, $v2_d)
```

    ┌ Warning: `Target(triple::String)` is deprecated, use `Target(; triple = triple)` instead.
    │   caller = ip:0x0
    └ @ Core :-1





    BenchmarkTools.Trial: 
      memory estimate:  3.28 KiB
      allocs estimate:  130
      --------------
      minimum time:     22.319 ms (0.00% GC)
      median time:      22.429 ms (0.00% GC)
      mean time:        22.722 ms (0.00% GC)
      maximum time:     26.938 ms (0.00% GC)
      --------------
      samples:          220
      evals/sample:     1



For CuSnpArray, the additional cost for mean imputation is negligible.


```julia
@benchmark CUDA.@sync LinearAlgebra.mul!($v1_d_, $EUR_100_cu_, $v2_d)
```




    BenchmarkTools.Trial: 
      memory estimate:  3.28 KiB
      allocs estimate:  130
      --------------
      minimum time:     22.243 ms (0.00% GC)
      median time:      22.431 ms (0.00% GC)
      mean time:        22.805 ms (0.00% GC)
      maximum time:     57.050 ms (0.00% GC)
      --------------
      samples:          220
      evals/sample:     1




```julia
EUR_100_mat_d = adapt(CuArray, EUR_100_mat);
```


```julia
@benchmark CUDA.@sync LinearAlgebra.mul!($v1_d, $EUR_100_mat_d, $v2_d)
```




    BenchmarkTools.Trial: 
      memory estimate:  2.58 KiB
      allocs estimate:  85
      --------------
      minimum time:     75.951 ms (0.00% GC)
      median time:      76.516 ms (0.00% GC)
      mean time:        77.909 ms (0.00% GC)
      maximum time:     82.096 ms (0.00% GC)
      --------------
      samples:          65
      evals/sample:     1



The speedup is obvious, CuSnpArrays is 30-50x faster than on CPU, and using CuSnpArray is both faster and memory-efficient compared to linear algebra with floating point matrix on GPU.


```julia
isapprox(v1_d, v1_d_)
```




    true



## $y = A^T x$


```julia
v1 = randn(size(EUR_100, 1))
v2 = randn(size(EUR_100, 2))
v1_d = adapt(CuArray{Float64}, v1)
v2_d = adapt(CuArray{Float64}, v2);
```


```julia
@benchmark LinearAlgebra.mul!($v2, transpose($EUR_100_sla), $v1)
```




    BenchmarkTools.Trial: 
      memory estimate:  16 bytes
      allocs estimate:  1
      --------------
      minimum time:     842.187 ms (0.00% GC)
      median time:      882.157 ms (0.00% GC)
      mean time:        913.111 ms (0.00% GC)
      maximum time:     1.049 s (0.00% GC)
      --------------
      samples:          6
      evals/sample:     1




```julia
@benchmark (LinearAlgebra.mul!($v2, transpose($EUR_100_bm), $v1))
```




    BenchmarkTools.Trial: 
      memory estimate:  16 bytes
      allocs estimate:  1
      --------------
      minimum time:     618.434 ms (0.00% GC)
      median time:      630.293 ms (0.00% GC)
      mean time:        689.063 ms (0.00% GC)
      maximum time:     881.782 ms (0.00% GC)
      --------------
      samples:          8
      evals/sample:     1




```julia
@benchmark LinearAlgebra.mul!($v2_d, transpose($EUR_100_cu), $v1_d)
```




    BenchmarkTools.Trial: 
      memory estimate:  3.08 KiB
      allocs estimate:  118
      --------------
      minimum time:     26.710 ms (0.00% GC)
      median time:      26.930 ms (0.00% GC)
      mean time:        27.350 ms (0.00% GC)
      maximum time:     32.306 ms (0.00% GC)
      --------------
      samples:          183
      evals/sample:     1




```julia
isapprox(collect(v2_d), v2)
```




    true




```julia
v1 = randn(size(EUR_101, 1))
v2 = randn(size(EUR_101, 2));
```


```julia
@benchmark LinearAlgebra.mul!($v2, transpose($EUR_101_sla), $v1)
```




    BenchmarkTools.Trial: 
      memory estimate:  128 bytes
      allocs estimate:  3
      --------------
      minimum time:     850.050 ms (0.00% GC)
      median time:      862.818 ms (0.00% GC)
      mean time:        953.739 ms (0.00% GC)
      maximum time:     1.391 s (0.00% GC)
      --------------
      samples:          6
      evals/sample:     1




```julia
@benchmark LinearAlgebra.mul!($v2, transpose($EUR_101_sla_), $v1)
```




    BenchmarkTools.Trial: 
      memory estimate:  128 bytes
      allocs estimate:  3
      --------------
      minimum time:     673.269 ms (0.00% GC)
      median time:      693.418 ms (0.00% GC)
      mean time:        700.729 ms (0.00% GC)
      maximum time:     784.130 ms (0.00% GC)
      --------------
      samples:          8
      evals/sample:     1




```julia
@benchmark (LinearAlgebra.mul!($v2, transpose($EUR_101_bm), $v1))
```




    BenchmarkTools.Trial: 
      memory estimate:  16 bytes
      allocs estimate:  1
      --------------
      minimum time:     620.830 ms (0.00% GC)
      median time:      637.698 ms (0.00% GC)
      mean time:        656.708 ms (0.00% GC)
      maximum time:     798.543 ms (0.00% GC)
      --------------
      samples:          8
      evals/sample:     1



BitMatrix is slightly faster in this direction.
