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
EUR_100_mat = convert(Matrix{Float64}, EUR_100, model=ADDITIVE_MODEL, center=true, scale=true);

EUR_101_bm = SnpBitMatrix{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_101_sla = SnpLinAlg{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false)
EUR_101_mat = convert(Matrix{Float64}, EUR_101, model=ADDITIVE_MODEL, center=true, scale=true);
```


```julia
ENV["JULIA_CUDA_USE_BINARYBUILDER"] = "false"
using CUDA
EUR_100_cu = CuSnpArray{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false);
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
      minimum time:     321.207 ms (0.00% GC)
      median time:      347.112 ms (0.00% GC)
      mean time:        375.497 ms (0.00% GC)
      maximum time:     769.390 ms (0.00% GC)
      --------------
      samples:          15
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
      minimum time:     1.745 s (0.00% GC)
      median time:      1.905 s (0.00% GC)
      mean time:        2.144 s (0.00% GC)
      maximum time:     2.782 s (0.00% GC)
      --------------
      samples:          3
      evals/sample:     1



Direct linear algebra on a SnpArray: 


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_100_sla, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  8.23 KiB
      allocs estimate:  158
      --------------
      minimum time:     1.014 s (0.00% GC)
      median time:      1.021 s (0.00% GC)
      mean time:        1.020 s (0.00% GC)
      maximum time:     1.027 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1



The below is the benchmark for SnpBitMatrix:


```julia
@benchmark (LinearAlgebra.mul!($v1, $EUR_100_bm, $v2))
```




    BenchmarkTools.Trial: 
      memory estimate:  0 bytes
      allocs estimate:  0
      --------------
      minimum time:     1.041 s (0.00% GC)
      median time:      1.056 s (0.00% GC)
      mean time:        1.058 s (0.00% GC)
      maximum time:     1.076 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1



Now let's try CUDA. The device is Nvidia Titan V.


```julia
using Adapt
```

Moving data to GPU: 


```julia
v1_d = adapt(CuArray{Float64}, v1)
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
      minimum time:     22.138 ms (0.00% GC)
      median time:      22.356 ms (0.00% GC)
      mean time:        22.352 ms (0.00% GC)
      maximum time:     23.825 ms (0.00% GC)
      --------------
      samples:          224
      evals/sample:     1



The speedup is obvious. Let's check correctness:


```julia
isapprox(collect(v1_d), v1)
```




    true




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
      minimum time:     78.002 ms (0.00% GC)
      median time:      80.142 ms (0.00% GC)
      mean time:        80.055 ms (0.00% GC)
      maximum time:     83.111 ms (0.00% GC)
      --------------
      samples:          63
      evals/sample:     1



Using CuSnpArray is both faster and memory-efficient compared to linear algebra with floating point matrix on GPU.

At first glance, the result from SnpBitMatrix might look similar to SnpLinAlg. However, SnpLinAlg is more stable in performance when the number of samples is not multiple of 4 or 8.


```julia
v1 = randn(size(EUR_101, 1))
v2 = randn(size(EUR_101, 2));
```


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_sla, $v2)
```




    BenchmarkTools.Trial: 
      memory estimate:  14.03 KiB
      allocs estimate:  264
      --------------
      minimum time:     1.022 s (0.00% GC)
      median time:      1.024 s (0.00% GC)
      mean time:        1.024 s (0.00% GC)
      maximum time:     1.028 s (0.00% GC)
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
      minimum time:     1.187 s (0.00% GC)
      median time:      1.188 s (0.00% GC)
      mean time:        1.191 s (0.00% GC)
      maximum time:     1.197 s (0.00% GC)
      --------------
      samples:          5
      evals/sample:     1



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
      minimum time:     662.131 ms (0.00% GC)
      median time:      678.294 ms (0.00% GC)
      mean time:        701.543 ms (0.00% GC)
      maximum time:     790.350 ms (0.00% GC)
      --------------
      samples:          8
      evals/sample:     1




```julia
@benchmark (LinearAlgebra.mul!($v2, transpose($EUR_100_bm), $v1))
```




    BenchmarkTools.Trial: 
      memory estimate:  16 bytes
      allocs estimate:  1
      --------------
      minimum time:     605.351 ms (0.00% GC)
      median time:      618.706 ms (0.00% GC)
      mean time:        620.553 ms (0.00% GC)
      maximum time:     658.692 ms (0.00% GC)
      --------------
      samples:          9
      evals/sample:     1




```julia
@benchmark LinearAlgebra.mul!($v2_d, transpose($EUR_100_cu), $v1_d)
```




    BenchmarkTools.Trial: 
      memory estimate:  3.08 KiB
      allocs estimate:  118
      --------------
      minimum time:     26.926 ms (0.00% GC)
      median time:      27.330 ms (0.00% GC)
      mean time:        27.494 ms (0.00% GC)
      maximum time:     31.492 ms (0.00% GC)
      --------------
      samples:          182
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
      minimum time:     667.369 ms (0.00% GC)
      median time:      681.264 ms (0.00% GC)
      mean time:        697.866 ms (0.00% GC)
      maximum time:     789.988 ms (0.00% GC)
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
      minimum time:     612.441 ms (0.00% GC)
      median time:      624.015 ms (0.00% GC)
      mean time:        641.973 ms (0.00% GC)
      maximum time:     770.845 ms (0.00% GC)
      --------------
      samples:          8
      evals/sample:     1



BitMatrix is slightly faster in this direction.
