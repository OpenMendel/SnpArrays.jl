# Linear Algebra Benchmarks

`SnpArrays.jl` supports three modes of matrix-vector multiplications:

`SnpArrays.jl` supports three modes of matrix-vector multiplications:

1. Direct operations on a plink-formatted `SnpArray`: `SnpLinAlg`
2. Operations on transformed `BitMatrix`es: `SnpBitMatrix` (support for this will be dropped in the near future)
3. Direct operations on a plink-formatted data on an Nvidia GPU: `CuSnpArray`.

`SnpLinAlg` also supports matrix-matrix multiplications.

- `SnpLinAlg` and `SnpBitMatrix` use Chris Elrod's [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl) internally. It is much faster on machines with AVX support.  
- `CuSnpArray` uses [CUDA.jl](https://juliagpu.gitlab.io/CUDA.jl/) internally.
On this page, we compare these three.
- `SnpLinAlg` supports multithreading. See [this page](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads-1) to learn how to use it.


```julia
versioninfo()
```

    Julia Version 1.6.2
    Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin19.6.0)
      CPU: Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-11.0.1 (ORCJIT, skylake-avx512)



```julia
using SnpArrays
using LinearAlgebra
using BenchmarkTools
```

    ┌ Info: Precompiling SnpArrays [4e780e97-f5bf-4111-9dc4-b70aaf691b06]
    └ @ Base loading.jl:1317
    WARNING: could not import DataFrames.DataFrame! into SnpArrays



```julia
const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"));
```

Let's try with EUR data repeated 100 and 101 times: 37900 by 54051 and 38279 by 54051, respectively.


```julia
EUR_10 = [EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR]
EUR_100 = [EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10; EUR_10];
EUR_101 = [EUR_100; EUR];
```

We create instances of SnpLinAlg, SnpBitmatrix and CuSnpArray:


```julia
EUR_100_bm = SnpBitMatrix{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_100_sla = SnpLinAlg{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_100_sla_ = SnpLinAlg{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false, impute=false);
EUR_100_mat = convert(Matrix{Float64}, EUR_100, model=ADDITIVE_MODEL, center=false, scale=false);

EUR_101_bm = SnpBitMatrix{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_101_sla = SnpLinAlg{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_101_sla_ = SnpLinAlg{Float64}(EUR_101; model=ADDITIVE_MODEL, center=false, scale=false, impute=false);
EUR_101_mat = convert(Matrix{Float64}, EUR_101, model=ADDITIVE_MODEL, center=false, scale=false);
```


```julia
using CUDA
EUR_100_cu = CuSnpArray{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false);
EUR_100_cu_ = CuSnpArray{Float64}(EUR_100; model=ADDITIVE_MODEL, center=false, scale=false, impute=false);
```

    ┌ Warning: The NVIDIA driver on this system only supports up to CUDA 10.2.0.
    │ For performance reasons, it is recommended to upgrade to a driver that supports CUDA 11.2 or higher.
    └ @ CUDA /home/xyz/.julia/packages/CUDA/CtvPY/src/initialization.jl:42


## $y = Ax$


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




    BenchmarkTools.Trial: 11 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m408.568 ms[22m[39m … [35m478.025 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m461.873 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m457.339 ms[22m[39m ± [32m 19.116 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m█[39m [39m [39m [39m█[39m [39m [39m [32m [39m[39m [34m█[39m[39m█[39m [39m [39m [39m█[39m█[39m [39m█[39m [39m [39m█[39m [39m [39m [39m█[39m [39m 
      [39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[32m▁[39m[39m▁[34m█[39m[39m█[39m▁[39m▁[39m▁[39m█[39m█[39m▁[39m█[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m█[39m [39m▁
      409 ms[90m           Histogram: frequency by time[39m          478 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



With single-threaded OpenBLAS: 


```julia
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($v1, $EUR_100_mat, $v2)
```




    BenchmarkTools.Trial: 3 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m1.957 s[22m[39m … [35m  2.089 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m1.986 s              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m2.011 s[22m[39m ± [32m69.424 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      1.96 s[90m         Histogram: frequency by time[39m         290 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



Direct linear algebra on a SnpArray, with mean imputation: 


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_100_sla, $v2)
```




    BenchmarkTools.Trial: 6 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m840.116 ms[22m[39m … [35m852.731 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m844.339 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m845.643 ms[22m[39m ± [32m  4.589 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m█[34m [39m[39m [39m [39m [39m█[39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m█[34m▁[39m[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      840 ms[90m           Histogram: frequency by time[39m          853 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m160 bytes[39m, allocs estimate[90m: [39m[33m1[39m.



With zero imputation:


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_100_sla_, $v2)
```




    BenchmarkTools.Trial: 9 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m590.436 ms[22m[39m … [35m600.609 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m594.370 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m594.687 ms[22m[39m ± [32m  2.947 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m▁[39m [39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m█[34m [39m[32m [39m[39m [39m [39m [39m▁[39m▁[39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[34m▁[39m[32m▁[39m[39m▁[39m▁[39m▁[39m█[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      590 ms[90m           Histogram: frequency by time[39m          601 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m160 bytes[39m, allocs estimate[90m: [39m[33m1[39m.



Indeed, we are paying some price for mean imputation.

The below is the benchmark for SnpBitMatrix (always zero-imputed):


```julia
@benchmark (LinearAlgebra.mul!($v1, $EUR_100_bm, $v2))
```




    BenchmarkTools.Trial: 2 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m4.945 s[22m[39m … [35m  4.981 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m4.963 s              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m4.963 s[22m[39m ± [32m25.520 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      4.94 s[90m         Histogram: frequency by time[39m        4.98 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



At first glance, the result from SnpBitMatrix might look better than SnpLinAlg. However, SnpLinAlg is more stable in performance when the number of samples is not multiple of 4 or 8.


```julia
v1 = randn(size(EUR_101, 1))
v2 = randn(size(EUR_101, 2));
```


```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_sla, $v2)
```




    BenchmarkTools.Trial: 6 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m858.307 ms[22m[39m … [35m895.131 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m867.094 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m869.931 ms[22m[39m ± [32m 13.289 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m█[39m [39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m [34m█[39m[39m [39m [39m [39m█[39m [39m [32m [39m[39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [39m█[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[34m█[39m[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[32m▁[39m[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      858 ms[90m           Histogram: frequency by time[39m          895 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m160 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_sla_, $v2)
```




    BenchmarkTools.Trial: 9 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m615.004 ms[22m[39m … [35m631.572 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m616.410 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m618.802 ms[22m[39m ± [32m  5.452 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m▁[39m▁[39m [39m [34m█[39m[39m█[39m [39m [39m [39m [39m [39m▁[39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [39m█[39m█[39m▁[39m▁[34m█[39m[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      615 ms[90m           Histogram: frequency by time[39m          632 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m160 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
@benchmark LinearAlgebra.mul!($v1, $EUR_101_bm, $v2)
```




    BenchmarkTools.Trial: 2 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m4.967 s[22m[39m … [35m  4.995 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m4.981 s              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m4.981 s[22m[39m ± [32m19.665 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      4.97 s[90m         Histogram: frequency by time[39m        4.99 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



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




    BenchmarkTools.Trial: 240 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m20.808 ms[22m[39m … [35m 23.129 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m20.834 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m20.857 ms[22m[39m ± [32m200.753 μs[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.08% ± 1.19%
    
      [39m▅[34m█[39m[32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m█[34m█[39m[32m█[39m[39m▄[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▄[39m [39m▅
      20.8 ms[90m       [39m[90mHistogram: [39m[90m[1mlog([22m[39m[90mfrequency[39m[90m[1m)[22m[39m[90m by time[39m      22.1 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m69.23 KiB[39m, allocs estimate[90m: [39m[33m4393[39m.



For CuSnpArray, the additional cost for mean imputation is negligible.


```julia
@benchmark CUDA.@sync LinearAlgebra.mul!($v1_d_, $EUR_100_cu_, $v2_d)
```




    BenchmarkTools.Trial: 239 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m20.813 ms[22m[39m … [35m 30.895 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m20.837 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m20.945 ms[22m[39m ± [32m915.783 μs[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.08% ± 1.27%
    
      [34m█[39m[32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [34m█[39m[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▄[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▄[39m [39m▅
      20.8 ms[90m       [39m[90mHistogram: [39m[90m[1mlog([22m[39m[90mfrequency[39m[90m[1m)[22m[39m[90m by time[39m      27.4 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m18.95 KiB[39m, allocs estimate[90m: [39m[33m1175[39m.




```julia
EUR_100_mat_d = adapt(CuArray, EUR_100_mat);
```


    Out of GPU memory trying to allocate 15.263 GiB
    Effective GPU memory usage: 11.12% (1.311 GiB/11.784 GiB)
    CUDA allocator usage: 980.317 MiB
    Memory pool usage: 980.317 MiB (980.317 MiB allocated, 0 bytes cached)


    

    Stacktrace:

      [1] #alloc#176

        @ ~/.julia/packages/CUDA/CtvPY/src/pool.jl:267 [inlined]

      [2] alloc

        @ ~/.julia/packages/CUDA/CtvPY/src/pool.jl:259 [inlined]

      [3] CuArray{Float64, 2}(#unused#::UndefInitializer, dims::Tuple{Int64, Int64})

        @ CUDA ~/.julia/packages/CUDA/CtvPY/src/array.jl:28

      [4] CuArray

        @ ~/.julia/packages/CUDA/CtvPY/src/array.jl:241 [inlined]

      [5] CuArray

        @ ~/.julia/packages/CUDA/CtvPY/src/array.jl:249 [inlined]

      [6] convert

        @ ~/.julia/packages/GPUArrays/Tebtl/src/host/construction.jl:4 [inlined]

      [7] adapt_storage(#unused#::Type{CuArray}, xs::Matrix{Float64})

        @ CUDA ~/.julia/packages/CUDA/CtvPY/src/array.jl:286

      [8] adapt_structure(to::Type, x::Matrix{Float64})

        @ Adapt ~/.julia/packages/Adapt/RGNRk/src/Adapt.jl:42

      [9] adapt(to::Type, x::Matrix{Float64})

        @ Adapt ~/.julia/packages/Adapt/RGNRk/src/Adapt.jl:40

     [10] top-level scope

        @ In[21]:1

     [11] eval

        @ ./boot.jl:360 [inlined]

     [12] include_string(mapexpr::typeof(REPL.softscope), mod::Module, code::String, filename::String)

        @ Base ./loading.jl:1116



```julia
@benchmark CUDA.@sync LinearAlgebra.mul!($v1_d, $EUR_100_mat_d, $v2_d)
```


    UndefVarError: EUR_100_mat_d not defined

    

    Stacktrace:

     [1] top-level scope

       @ ~/.julia/packages/BenchmarkTools/tGTCy/src/execution.jl:440

     [2] eval

       @ ./boot.jl:360 [inlined]

     [3] include_string(mapexpr::typeof(REPL.softscope), mod::Module, code::String, filename::String)

       @ Base ./loading.jl:1116


The speedup is obvious, CuSnpArrays is 30-50x faster than on CPU, and using CuSnpArray is both faster and memory-efficient compared to linear algebra with floating point matrix on GPU.


```julia
isapprox(v1_d, v1_d_)
```




    true



## $y = A^T x$


```julia
v1 = randn(size(EUR_100, 1))
v2 = randn(size(EUR_100, 2))
v2_ = randn(size(EUR_100, 2))
v1_d = adapt(CuArray{Float64}, v1)
v2_d = adapt(CuArray{Float64}, v2);
```


```julia
@benchmark LinearAlgebra.mul!($v2, transpose($EUR_100_sla), $v1)
```




    BenchmarkTools.Trial: 6 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m851.938 ms[22m[39m … [35m875.831 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m864.144 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m864.772 ms[22m[39m ± [32m  8.569 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m [39m█[34m█[39m[39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[34m█[39m[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      852 ms[90m           Histogram: frequency by time[39m          876 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m304 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
@benchmark (LinearAlgebra.mul!($v2, transpose($EUR_100_bm), $v1))
```




    BenchmarkTools.Trial: 1 sample with 1 evaluation.
     Single result which took [34m5.943 s[39m (0.00% GC) to evaluate,
     with a memory estimate of [33m0 bytes[39m, over [33m0[39m allocations.




```julia
@benchmark LinearAlgebra.mul!($v2_d, transpose($EUR_100_cu), $v1_d)
```




    BenchmarkTools.Trial: 278 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m17.848 ms[22m[39m … [35m52.127 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 65.80%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m17.878 ms              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m18.038 ms[22m[39m ± [32m 2.084 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.68% ±  3.95%
    
      [39m█[34m▅[39m[39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m█[34m█[39m[39m▆[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▄[39m [39m▅
      17.8 ms[90m      [39m[90mHistogram: [39m[90m[1mlog([22m[39m[90mfrequency[39m[90m[1m)[22m[39m[90m by time[39m        20 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m142.14 KiB[39m, allocs estimate[90m: [39m[33m9059[39m.




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




    BenchmarkTools.Trial: 5 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m868.448 ms[22m[39m … [35m   1.327 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m895.181 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m   1.033 s[22m[39m ± [32m217.498 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [34m█[39m[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      868 ms[90m           Histogram: frequency by time[39m          1.33 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m304 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
@benchmark LinearAlgebra.mul!($v2, transpose($EUR_101_sla_), $v1)
```




    BenchmarkTools.Trial: 6 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m927.944 ms[22m[39m … [35m  1.028 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m934.817 ms              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m958.481 ms[22m[39m ± [32m42.785 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m▁[39m█[34m [39m[39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [39m█[39m█[34m▁[39m[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      928 ms[90m          Histogram: frequency by time[39m           130 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m304 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
@benchmark (LinearAlgebra.mul!($v2, transpose($EUR_101_bm), $v1))
```




    BenchmarkTools.Trial: 1 sample with 1 evaluation.
     Single result which took [34m5.813 s[39m (0.00% GC) to evaluate,
     with a memory estimate of [33m0 bytes[39m, over [33m0[39m allocations.



BitMatrix is slightly faster in this direction.

## $Y = AX$

Now for matrix-matrix multiplications. If we want to center/scale the SnpArray, we have
```math
\begin{aligned}
    Y_{ij} = \sum_{k} \frac{A_{ik} - \mu_k}{\sigma_k}X_{kj}
\end{aligned}
```

where $\mu_k$ and $\sigma_k$ is the mean and standard deviation of the $k$th SNP. Centering and scaling is performed on-the-fly. First check correctness


```julia
EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"));
EUR_10 = [EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR]

m = size(EUR_10, 1)
n = size(EUR_10, 2)
p = 2

A = SnpLinAlg{Float64}(EUR_10; model=ADDITIVE_MODEL, impute=true, center=true, scale=true);
X = rand(n, p)
Y = zeros(m, p)
SnpArrays.mul!(Y, A, X)
Afloat = convert(Matrix{Float64}, EUR_10, impute=true, center=true, scale=true)
Ytrue = Afloat * X
all(Y .≈ Ytrue)
```

    WARNING: redefinition of constant EUR. This may fail, cause incorrect answers, or produce other errors.





    true



Now lets check out timings. If $B$ is a "tall and thin" matrix, then `SnpLinAlg` remains competitive, often superior, to BLAS. 


```julia
# SnpLinAlg-matrix
@benchmark LinearAlgebra.mul!($Y, $A, $X)
```




    BenchmarkTools.Trial: 46 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m105.415 ms[22m[39m … [35m132.895 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m106.670 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m108.777 ms[22m[39m ± [32m  6.331 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m [39m▄[39m█[34m▅[39m[39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m▇[39m█[39m█[34m█[39m[39m▃[39m▃[39m▁[32m▁[39m[39m▃[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▃[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▃[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▄[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▃[39m▁[39m▃[39m [39m▁
      105 ms[90m           Histogram: frequency by time[39m          133 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m112 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
# BLAS with 1 threaed 
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($Y, $Afloat, $X)
```




    BenchmarkTools.Trial: 8 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m589.561 ms[22m[39m … [35m885.109 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m601.005 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m635.873 ms[22m[39m ± [32m100.867 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m [39m [34m█[39m[39m▃[39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m▇[39m▇[34m█[39m[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▇[39m [39m▁
      590 ms[90m           Histogram: frequency by time[39m          885 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.




```julia
# BLAS with 8 threaed 
BLAS.set_num_threads(8)
@benchmark LinearAlgebra.mul!($Y, $Afloat, $X)
```




    BenchmarkTools.Trial: 40 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m109.089 ms[22m[39m … [35m609.969 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m111.358 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m127.126 ms[22m[39m ± [32m 78.665 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [34m█[39m[39m▃[32m▃[39m[39m▃[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▂[39m [39m▁
      109 ms[90m           Histogram: frequency by time[39m          610 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



But if $B$ is a large matrix too, single threaded BLAS is ~6 times faster and >10x faster for 8 thread BLAS. 


```julia
# SnpLinAlg-matrix
p = 100
X = rand(n, p)
Y = zeros(m, p)
@benchmark LinearAlgebra.mul!($Y, $A, $X)
```




    BenchmarkTools.Trial: 2 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m3.599 s[22m[39m … [35m 3.602 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m3.600 s             [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m3.600 s[22m[39m ± [32m1.921 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      3.6 s[90m         Histogram: frequency by time[39m         3.6 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m112 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
# BLAS with 1 threaed 
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($Y, $Afloat, $X)
```




    BenchmarkTools.Trial: 3 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m2.192 s[22m[39m … [35m   2.759 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m2.192 s               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m2.381 s[22m[39m ± [32m326.877 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      2.19 s[90m         Histogram: frequency by time[39m         2.76 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.




```julia
# BLAS with 8 threaed 
BLAS.set_num_threads(8)
@benchmark LinearAlgebra.mul!($Y, $Afloat, $X)
```




    BenchmarkTools.Trial: 15 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m323.403 ms[22m[39m … [35m366.700 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m330.332 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m334.467 ms[22m[39m ± [32m 11.791 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m▁[39m█[39m▁[39m [39m [39m█[34m▁[39m[39m [39m [39m▁[39m [39m [39m [39m [39m [39m [32m▁[39m[39m▁[39m [39m [39m [39m▁[39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m▁[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m▁[39m [39m 
      [39m█[39m█[39m█[39m▁[39m▁[39m█[34m█[39m[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m█[39m[39m█[39m▁[39m▁[39m▁[39m█[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m█[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      323 ms[90m           Histogram: frequency by time[39m          367 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



## $Y = A^tX$

If we want to center/scale the SnpArray, we have

```math
\begin{aligned}
    Y_{ij} = \sum_{k} \left(\frac{A_{ik} - \mu_i}{\sigma_i}\right)X_{kj}
\end{aligned}
```

where $\mu_i$ and $\sigma_i$ is the mean and standard deviation of the $i$th SNP. Similar to before, lets first check correctness.


```julia
EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"));
EUR_10 = [EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR; EUR]

m = size(EUR_10, 1)
n = size(EUR_10, 2)
p = 2

A = SnpLinAlg{Float64}(EUR_10; model=ADDITIVE_MODEL, impute=true, center=true, scale=true);
X = rand(m, p)
Y = zeros(n, p)
SnpArrays.mul!(Y, Transpose(A), X)
Afloat = convert(Matrix{Float64}, EUR_10, impute=true, center=true, scale=true)
Ytrue = Afloat' * X
all(Y .≈ Ytrue)
```

    WARNING: redefinition of constant EUR. This may fail, cause incorrect answers, or produce other errors.





    true



Now lets check out timings. If $B$ is a "tall and thin" matrix, then `SnpLinAlg` remains competitive, often superior, to BLAS. 


```julia
# SnpLinAlg-matrix
@benchmark LinearAlgebra.mul!($Y, $(Transpose(A)), $X)
```




    BenchmarkTools.Trial: 25 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m193.170 ms[22m[39m … [35m217.577 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m199.180 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m200.356 ms[22m[39m ± [32m  5.605 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m [39m [39m [39m▃[39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m▃[34m [39m[39m [39m [39m▃[32m▃[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m▇[39m▁[39m▁[39m█[39m▁[39m▇[39m▇[39m▇[39m▁[39m▇[39m█[39m▁[39m▁[39m█[34m▁[39m[39m▇[39m▁[39m█[32m█[39m[39m▇[39m▁[39m▇[39m▇[39m▁[39m▁[39m▁[39m▁[39m▇[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▇[39m▁[39m▁[39m▁[39m▇[39m▁[39m▁[39m▁[39m▇[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▇[39m [39m▁
      193 ms[90m           Histogram: frequency by time[39m          218 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m304 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
# BLAS with 1 threaed 
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($Y, $(Transpose(Afloat)), $X)
```




    BenchmarkTools.Trial: 13 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m402.128 ms[22m[39m … [35m452.548 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m404.967 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m408.729 ms[22m[39m ± [32m 13.282 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m [39m█[39m [34m█[39m[39m▃[39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m▇[39m█[39m▁[34m█[39m[39m█[39m▇[39m▇[39m▇[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▇[39m [39m▁
      402 ms[90m           Histogram: frequency by time[39m          453 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.




```julia
# BLAS with 8 threaed 
BLAS.set_num_threads(8)
@benchmark LinearAlgebra.mul!($Y, $(Transpose(Afloat)), $X)
```




    BenchmarkTools.Trial: 74 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m65.019 ms[22m[39m … [35m96.348 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m65.594 ms              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m67.762 ms[22m[39m ± [32m 5.614 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m▅[39m█[34m▂[39m[39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m█[39m█[34m█[39m[39m▁[39m▅[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▅[39m▅[39m▆[39m▅[39m▅[39m▅[39m▁[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m▅[39m▁[39m▁[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m▅[39m [39m▁
      65 ms[90m        [39m[90mHistogram: [39m[90m[1mlog([22m[39m[90mfrequency[39m[90m[1m)[22m[39m[90m by time[39m      86.5 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



But if $B$ is a large matrix too, single threaded BLAS is ~6 times faster and >10x faster for 8 thread BLAS. 


```julia
# SnpLinAlg-matrix
p = 100
X = rand(m, p)
Y = zeros(n, p)
@benchmark LinearAlgebra.mul!($Y, $(Transpose(A)), $X)
```




    BenchmarkTools.Trial: 2 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m4.218 s[22m[39m … [35m   4.367 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m4.293 s               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m4.293 s[22m[39m ± [32m105.845 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      4.22 s[90m         Histogram: frequency by time[39m         4.37 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m304 bytes[39m, allocs estimate[90m: [39m[33m1[39m.




```julia
# BLAS with 1 threaed 
BLAS.set_num_threads(1)
@benchmark LinearAlgebra.mul!($Y, $(Transpose(Afloat)), $X)
```




    BenchmarkTools.Trial: 3 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m1.845 s[22m[39m … [35m  1.870 s[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m1.855 s              [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m1.857 s[22m[39m ± [32m12.738 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [34m█[39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m█[39m [39m 
      [34m█[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m▁[39m▁[39m▁[32m▁[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m█[39m [39m▁
      1.84 s[90m         Histogram: frequency by time[39m        1.87 s [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.




```julia
# BLAS with 8 threaed 
BLAS.set_num_threads(8)
@benchmark LinearAlgebra.mul!($Y, $(Transpose(Afloat)), $X)
```




    BenchmarkTools.Trial: 17 samples with 1 evaluation.
     Range [90m([39m[36m[1mmin[22m[39m … [35mmax[39m[90m):  [39m[36m[1m276.946 ms[22m[39m … [35m386.706 ms[39m  [90m┊[39m GC [90m([39mmin … max[90m): [39m0.00% … 0.00%
     Time  [90m([39m[34m[1mmedian[22m[39m[90m):     [39m[34m[1m285.168 ms               [22m[39m[90m┊[39m GC [90m([39mmedian[90m):    [39m0.00%
     Time  [90m([39m[32m[1mmean[22m[39m ± [32mσ[39m[90m):   [39m[32m[1m302.769 ms[22m[39m ± [32m 37.051 ms[39m  [90m┊[39m GC [90m([39mmean ± σ[90m):  [39m0.00% ± 0.00%
    
      [39m█[39m [39m [34m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [32m [39m[39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m [39m 
      [39m█[39m▅[39m█[34m▁[39m[39m▅[39m▅[39m▁[39m▁[39m▁[39m▁[39m▅[39m▁[39m▁[39m▅[32m▅[39m[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m▁[39m▅[39m▁[39m▁[39m▁[39m▁[39m▁[39m▁[39m▅[39m [39m▁
      277 ms[90m           Histogram: frequency by time[39m          387 ms [0m[1m<[22m
    
     Memory estimate[90m: [39m[33m0 bytes[39m, allocs estimate[90m: [39m[33m0[39m.



## Conclusion

+ `SnpLinAlg` (for CPU)
    - achieves up to 32x memory savings compared to double-precision matrices
    - is usually faster than single threaded BLAS for matrix-vector multiply.
    - is competitive with single threaded BLAS for matrix-matrix multiply if $B$ is "tall and thin"
+ `CuSnpArray` supports GPU matrix-vector operations that is 30-50x faster than multithreaded BLAS. 
+ Other linear algebra operations (e.g. $v*A$ and $qr(A)$...etc) will be *much* slower and are not guaranteed to work. 
