# SnpArrays.jl

Data from [*genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) are often saved as a [**PLINK binary biallelic genotype table**](https://www.cog-genomics.org/plink2/formats#bed) or `.bed` file. To be useful, such files should be accompanied by a `.fam` file, containing metadata on the rows of the table, and a `.bim` file,
containing metadata on the columns. The `.fam` and `.bim` files are in tab-separated format.

The table contains the observed allelic type at `n` [*single nucleotide polymorphism*](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (SNP) positions for `m` individuals. A SNP corresponds to a nucleotide position on the genome where some degree of variation has been observed in a population, with each individual have one of two possible *alleles* at that position on each of a pair of chromosomes. Three possible genotypes and corresponding coding are

| Genotype | Plink/SnpArray |  
|:---:|:---:|  
| A1,A1 | 0x00 |  
| missing | 0x01 |
| A1,A2 | 0x02 |  
| A2,A2 | 0x03 |  

## Installation

This package requires Julia v1.4 or later, which can be obtained from
<https://julialang.org/downloads/> or by building Julia from the sources in the
<https://github.com/JuliaLang/julia> repository.

The package can be installed by running the following code:
```julia
using Pkg
pkg"add SnpArrays"
```
For running the examples below, the following are also necessary. 
```julia
pkg"add BenchmarkTools DelimitedFiles Glob"
pkg"add https://github.com/OpenMendel/ADMIXTURE.jl"
```

For optional use on a CUDA-enabled GPU, the following is also needed. 
```julia
pkg"add Adapt CUDA"
```


```julia
versioninfo()
```

    Julia Version 1.6.2
    Commit 1b93d53fc4 (2021-07-14 15:36 UTC)
    Platform Info:
      OS: Linux (x86_64-pc-linux-gnu)
      CPU: Intel(R) Xeon(R) Silver 4114 CPU @ 2.20GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-11.0.1 (ORCJIT, skylake-avx512)



```julia
# for use in this tutorial
using SnpArrays, ADMIXTURE, BenchmarkTools, DelimitedFiles, Glob
Sys.islinux() && (using CUDA);
```

## Example data

There are two example data sets attached to this package. They are availabe in the `data` folder of the package.


```julia
datapath = normpath(SnpArrays.datadir())
```




    "/home/xyz/.julia/dev/SnpArrays/data"




```julia
readdir(glob"mouse.*", datapath)
```




    3-element Vector{String}:
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.bed"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.bim"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.fam"



Data set `EUR_subset` contains no missing genotypes. It is located at


```julia
readdir(glob"EUR_subset.*", datapath)
```




    3-element Vector{String}:
     "/home/xyz/.julia/dev/SnpArrays/data/EUR_subset.bed"
     "/home/xyz/.julia/dev/SnpArrays/data/EUR_subset.bim"
     "/home/xyz/.julia/dev/SnpArrays/data/EUR_subset.fam"



Data from recent studies, which have samples from tens of thousands of individuals at over a million SNP positions, would be in the tens or even hundreds of Gb range.

## SnpArray

`SnpArray` is the fundamental type for dealing with genotype data in Plink bed file. Each row of `SnpArray` is a sample and each column a SNP.

### Constructor

There are various ways to initialize a SnpArray.

#### Intitialize from Plink file set

SnpArray can be initialized from the Plink bed file. The corresponding `.fam` needs to be present, which is used to determine the number of individuals.


```julia
const mouse = SnpArray(SnpArrays.datadir("mouse.bed"))
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                    
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



The virtual size of the GWAS data is 1940 observations at each of 10150 SNP positions.


```julia
size(mouse)
```




    (1940, 10150)



Because the file is memory-mapped opening the file and accessing the data is fast, even for very large .bed files.


```julia
@btime(SnpArray(SnpArrays.datadir("mouse.bed")));
```

      67.366 μs (57 allocations: 389.66 KiB)


By default, the memory-mapped file is read only, changing entries is not allowed.


```julia
mouse[1, 1] = 0x00
```


    ReadOnlyMemoryError()

    

    Stacktrace:

     [1] setindex!

       @ ./array.jl:841 [inlined]

     [2] setindex!(s::SnpArray, x::UInt8, i::Int64, j::Int64)

       @ SnpArrays ~/.julia/dev/SnpArrays/src/snparray.jl:131

     [3] top-level scope

       @ In[9]:1

     [4] eval

       @ ./boot.jl:360 [inlined]

     [5] include_string(mapexpr::typeof(REPL.softscope), mod::Module, code::String, filename::String)

       @ Base ./loading.jl:1116


To possibly change genoytpes in a bed file, open with write permission
```julia
mouse = SnpArray(SnpArrays.datadir("mouse.bed"), "w")
```

#### Initialize from only bed file

If only the bed file is present, user is required to supply the number of individuals in the second argument.


```julia
SnpArray(SnpArrays.datadir("mouse.bed"), 1940)
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                    
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



#### Initialize from compressed Plink files

SnpArray can be initialized from Plink files in compressed formats: `gz`, `zlib`, `zz`, `xz`, `zst`, or `bz2`. For a complete list type
```julia
SnpArrays.ALLOWED_FORMAT
```
If you want to support a new compressed format, file an issue.

Let us first compress the mouse data in gz format. We see gz format takes less than 1/3 storage of original Plink files.


```julia
compress_plink(SnpArrays.datadir("mouse"), "gz")
readdir(glob"mouse.*.gz", datapath)
```




    3-element Vector{String}:
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.bed.gz"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.bim.gz"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.fam.gz"



To initialize SnpArray from gzipped Plink file, simply used the bed file with name ending with `.bed.gz`:


```julia
# requires corresponding `.fam.gz` file
SnpArray(SnpArrays.datadir("mouse.bed.gz"))
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                    
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03



or


```julia
# does not require corresponding `.fam.gz` file
SnpArray(SnpArrays.datadir("mouse.bed.gz"), 1940)
```




    1940×10150 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x00  0x00  0x00  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x00  0x00  0x00  0x00  0x00  0x00
        ⋮                             ⋮  ⋱           ⋮                    
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x01  0x01  0x01  0x01  0x01  0x01
     0x00  0x00  0x00  0x00  0x03  0x00     0x03  0x03  0x03  0x03  0x03  0x03




```julia
# clean up
rm(SnpArrays.datadir("mouse.bed.gz"), force=true)
rm(SnpArrays.datadir("mouse.fam.gz"), force=true)
rm(SnpArrays.datadir("mouse.bim.gz"), force=true)
```

#### Initialize and create bed file

Initialize 5 rows and 3 columns with all (A1, A1) genotype (0x00) and memory-map to a bed file `tmp.bed`


```julia
tmpbf = SnpArray("tmp.bed", 5, 3)
```




    5×3 SnpArray:
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00



Change entries


```julia
tmpbf[1:2, 1:2] .= 0x03
tmpbf
```




    5×3 SnpArray:
     0x03  0x03  0x00
     0x03  0x03  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00
     0x00  0x00  0x00




```julia
fill!(tmpbf, 0x02)
tmpbf
```




    5×3 SnpArray:
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02
     0x02  0x02  0x02




```julia
# clean up
rm("tmp.bed", force=true)
```

Initialize 5 rows and 3 columns with undefined genotypes without memory-mapping to any file


```julia
tmpbf = SnpArray(undef, 5, 3)
```




    5×3 SnpArray:
     0x00  0x01  0x03
     0x00  0x02  0x01
     0x00  0x03  0x02
     0x00  0x03  0x00
     0x00  0x02  0x03



Create a bed file corresponding to an existing SnpArray and memory-map it.


```julia
tmpbf = SnpArray("tmp.bed", tmpbf)
```




    5×3 SnpArray:
     0x00  0x01  0x03
     0x00  0x02  0x01
     0x00  0x03  0x02
     0x00  0x03  0x00
     0x00  0x02  0x03




```julia
tmpbf[1, 1] = 0x02
tmpbf
```




    5×3 SnpArray:
     0x02  0x01  0x03
     0x00  0x02  0x01
     0x00  0x03  0x02
     0x00  0x03  0x00
     0x00  0x02  0x03




```julia
# clean up
rm("tmp.bed", force=true)
```

### `convert` and `copyto!`

Most common usage of SnpArray is to convert genotypes to numeric values for statistical analysis. Conversion rule depends on genetic models (additive, dominant, or recessive), centering, scaling, or imputation.

#### `convert`

`convert` function has 4 keyword arguments: `model`, `center`, `scale`, and `impute`.

`model` keyword specifies the SNP model for conversion. By default `convert` function translates genotypes according to the *additive* SNP model, which essentially counts the number of **A2** allele (0, 1 or 2) per genotype. Other SNP models are *dominant* and *recessive*, both in terms of the **A2** allele.

| Genotype | `SnpArray` | `model=ADDITIVE_MODEL` | `model=DOMINANT_MODEL` | `model=RECESSIVE_MODEL` |    
|:---:|:---:|:---:|:---:|:---:|  
| A1,A1 | 0x00 | 0 | 0 | 0 |  
| missing | 0x01 | NaN | NaN | NaN |
| A1,A2 | 0x02 | 1 | 1 | 0 |  
| A2,A2 | 0x03 | 2 | 1 | 1 |  

`center=true` tells `convert` to center each column by its mean. Default is `false`.

`scale=true` tells `convert` to scale each column by its standard deviation. Default is `false`.

`impute=true` tells `convert` to impute missing genotypes (0x01) by column mean. Default is `false`.

Convert whole SnpArray to a Float64 matrix using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`)


```julia
convert(Matrix{Float64}, mouse)
```




    1940×10150 Matrix{Float64}:
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       1.0    1.0    1.0    1.0    1.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       1.0    1.0    1.0    1.0    1.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0  …    0.0    0.0    0.0    0.0    0.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       0.0    0.0    0.0    0.0    0.0
     ⋮                        ⋮              ⋱    ⋮                         
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0       2.0    2.0    2.0    2.0    2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  1.0  …    2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0  2.0       2.0    2.0    2.0    2.0    2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  2.0     NaN    NaN    NaN    NaN    NaN
     0.0  0.0  0.0  0.0  2.0  0.0  2.0  0.0       2.0    2.0    2.0    2.0    2.0



!!! note  

    When `convert` or `copyto!` a slice or subarray of SnpArray, using `view`, `@view` or `views` is necessary for both correctness and efficiency. Without view, it's simply converting the UInt8 coding in original bed file.
    

Convert a column to Float64 vector using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`).


```julia
# convert(Vector{Float64}, view(mouse, :, 1)) # alternative syntax
# @views convert(Vector{Float64}, mouse[:, 1]) # alternative syntax
convert(Vector{Float64}, @view(mouse[:, 1]))
```




    1940-element Vector{Float64}:
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     2.0
     2.0
     1.0
     2.0
     ⋮
     2.0
     2.0
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     0.0



Convert a subarray of SnpArray to Float64 matrix using defaults (`model=ADDITIVE_MODEL`, `center=false`, `scale=false`, `impute=false`).


```julia
convert(Matrix{Float64}, @view(mouse[1:2:10, 1:2:10]))
```




    5×5 Matrix{Float64}:
     1.0  1.0  2.0  2.0  1.0
     2.0  2.0  2.0  2.0  2.0
     2.0  2.0  2.0  2.0  2.0
     1.0  1.0  2.0  2.0  1.0
     1.0  2.0  1.0  1.0  1.0



Different SNP models (`ADDITIVE_MODEL` vs `DOMINANT_MODEL` vs `RECESSIVE_MODEL`)


```julia
@views [convert(Vector{Float64}, mouse[:, 1], model=ADDITIVE_MODEL) convert(Vector{Float64}, mouse[:, 1], model=DOMINANT_MODEL) convert(Vector{Float64}, mouse[:, 1], model=RECESSIVE_MODEL)]
```




    1940×3 Matrix{Float64}:
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     ⋮         
     2.0  1.0  1.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     2.0  1.0  1.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     1.0  1.0  0.0
     0.0  0.0  0.0



Center and scale (last column) while `convert`


```julia
convert(Vector{Float64}, @view(mouse[:, end]), center=true, scale=true)
```




    1940-element Vector{Float64}:
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
      -1.8819155626127624
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
      -1.8819155626127624
      -4.2359771983402785
       0.4721460731147541
      -4.2359771983402785
       ⋮
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
       0.4721460731147541
     NaN
       0.4721460731147541



Center, scale, and impute (last column) while `convert`


```julia
convert(Vector{Float64}, @view(mouse[:, end]), center=true, scale=true, impute=true)
```




    1940-element Vector{Float64}:
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
     -1.8819155626127624
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
     -1.8819155626127624
     -4.2359771983402785
      0.4721460731147541
     -4.2359771983402785
      ⋮
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.4721460731147541
      0.0
      0.4721460731147541



#### `copyto!`

`copyto!` is the in-place version of `convert`. It takes the same keyword arguments (`model`, `center`, `scale`, `impute`) as `convert`.

Copy a column to a Float64 vector using defaults (`model=:additive`, `center=false`, `scale=false`, `impute=false`).


```julia
v = zeros(size(mouse, 1))
copyto!(v, @view(mouse[:, 1]))
```




    1940-element Vector{Float64}:
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     2.0
     2.0
     1.0
     2.0
     ⋮
     2.0
     2.0
     1.0
     1.0
     2.0
     1.0
     2.0
     1.0
     1.0
     1.0
     1.0
     0.0




```julia
@btime(copyto!($v, $@view(mouse[:, 1])));
```

      5.623 μs (0 allocations: 0 bytes)


Copy columns using defaults


```julia
v2 = zeros(size(mouse, 1), 2)
copyto!(v2, @view(mouse[:, 1:2]))
```




    1940×2 Matrix{Float64}:
     1.0  1.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     2.0  2.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     ⋮    
     2.0  2.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     2.0  2.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     1.0  1.0
     0.0  0.0




```julia
# roughly double the cost of copying 1 column
@btime(copyto!($v2, $@view(mouse[:, 1:2])));
```

      10.333 μs (0 allocations: 0 bytes)


Center and scale


```julia
copyto!(v, @view(mouse[:, 1]), center=true, scale=true)
```




    1940-element Vector{Float64}:
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626
     -0.16084075452851265
      1.2624897581484626
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626
      1.2624897581484626
     -0.16084075452851265
      1.2624897581484626
      ⋮
      1.2624897581484626
      1.2624897581484626
     -0.16084075452851265
     -0.16084075452851265
      1.2624897581484626
     -0.16084075452851265
      1.2624897581484626
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -0.16084075452851265
     -1.584171267205488




```julia
# more cost becoz of extra pass for center, scale, and/or impute
@btime(copyto!($v, $(@view(mouse[:, 1])), center=true, scale=true));
```

      6.381 μs (0 allocations: 0 bytes)


Looping over all columns


```julia
v = Vector{Float64}(undef, size(mouse, 1))
function loop_test(v, s)
    for j in 1:size(s, 2)
        copyto!(v, @view(s[:, j]))
    end
end
@btime(loop_test($v, $mouse))
```

      69.939 ms (0 allocations: 0 bytes)


Copy whole SnpArray


```julia
M = similar(mouse, Float64)
@btime(copyto!($M, $mouse));
```

      77.183 ms (0 allocations: 0 bytes)


#### Impute missing genotypes using ADMIXTURE estimates

`convert` and `copyto!` can perform more fine-tuned imputation using the ancestry estimates from the [ADMIXTURE](https://github.com/OpenMendel/ADMIXTURE.jl) software.

Step 1: Calculate the ancestry estimate and allele frequencies using ADMIXTURE.jl. Here we assume $K=3$ populations.


```julia
# install ADMIXTURE package first 
using ADMIXTURE
if isfile("mouse.3.P") && isfile("mouse.3.Q")
    P = readdlm("mouse.3.P", ' ', Float64) 
    Q = readdlm("mouse.3.Q", ' ', Float64)
else
    # run ADMIXTURE using 4 threads
    P, Q = admixture(SnpArrays.datadir("mouse.bed"), 3, j=4)
end;
```

    ****                   ADMIXTURE Version 1.3.0                  ****
    ****                    Copyright 2008-2015                     ****
    ****           David Alexander, Suyash Shringarpure,            ****
    ****                John  Novembre, Ken Lange                   ****
    ****                                                            ****
    ****                 Please cite our paper!                     ****
    ****   Information at www.genetics.ucla.edu/software/admixture  ****
    
    Parallel execution requested.  Will use 4 threads.
    Random seed: 43
    Point estimation method: Block relaxation algorithm
    Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
    Point estimation will terminate when objective function delta < 0.0001
    Estimation of standard errors disabled; will compute point estimates only.


    ┌ Info: ADMIXTURE command:
    │ `/home/xyz/.julia/artifacts/316b9c66aef8f67001d54aa86a244d1e769c1e1a/dist/admixture_linux-1.3.0/admixture /home/xyz/.julia/dev/SnpArrays/data/mouse.bed 3 -j4`
    └ @ ADMIXTURE /home/xyz/.julia/packages/ADMIXTURE/TuI9H/src/ADMIXTURE.jl:59
    ┌ Info: Output directory: /home/xyz
    └ @ ADMIXTURE /home/xyz/.julia/packages/ADMIXTURE/TuI9H/src/ADMIXTURE.jl:60


    Size of G: 1940x10150
    Performing five EM steps to prime main algorithm
    1 (EM) 	Elapsed: 1.434	Loglikelihood: -2.27484e+07	(delta): 8.92872e+06
    2 (EM) 	Elapsed: 1.14	Loglikelihood: -2.21886e+07	(delta): 559814
    3 (EM) 	Elapsed: 1.111	Loglikelihood: -2.20025e+07	(delta): 186060
    4 (EM) 	Elapsed: 1.133	Loglikelihood: -2.1896e+07	(delta): 106495
    5 (EM) 	Elapsed: 1.139	Loglikelihood: -2.18274e+07	(delta): 68590.1
    Initial loglikelihood: -2.18274e+07
    Starting main algorithm
    1 (QN/Block) 	Elapsed: 4.623	Loglikelihood: -2.12515e+07	(delta): 575921
    2 (QN/Block) 	Elapsed: 4.425	Loglikelihood: -2.10686e+07	(delta): 182932
    3 (QN/Block) 	Elapsed: 4.538	Loglikelihood: -2.09068e+07	(delta): 161743
    4 (QN/Block) 	Elapsed: 4.652	Loglikelihood: -2.07604e+07	(delta): 146489
    5 (QN/Block) 	Elapsed: 4.551	Loglikelihood: -2.07231e+07	(delta): 37298.4
    6 (QN/Block) 	Elapsed: 4.91	Loglikelihood: -2.07134e+07	(delta): 9625.64
    7 (QN/Block) 	Elapsed: 4.717	Loglikelihood: -2.07086e+07	(delta): 4869.1
    8 (QN/Block) 	Elapsed: 4.58	Loglikelihood: -2.07075e+07	(delta): 1085.47
    9 (QN/Block) 	Elapsed: 4.583	Loglikelihood: -2.07073e+07	(delta): 211.994
    10 (QN/Block) 	Elapsed: 4.59	Loglikelihood: -2.07072e+07	(delta): 39.2162
    11 (QN/Block) 	Elapsed: 4.585	Loglikelihood: -2.07072e+07	(delta): 8.32494
    12 (QN/Block) 	Elapsed: 4.586	Loglikelihood: -2.07072e+07	(delta): 1.5465
    13 (QN/Block) 	Elapsed: 4.538	Loglikelihood: -2.07072e+07	(delta): 0.168445
    14 (QN/Block) 	Elapsed: 4.576	Loglikelihood: -2.07072e+07	(delta): 0.0237868
    15 (QN/Block) 	Elapsed: 4.599	Loglikelihood: -2.07072e+07	(delta): 0.00178787
    16 (QN/Block) 	Elapsed: 4.622	Loglikelihood: -2.07072e+07	(delta): 0.000387926
    17 (QN/Block) 	Elapsed: 4.616	Loglikelihood: -2.07072e+07	(delta): 8.44523e-06
    Summary: 
    Converged in 17 iterations (85.726 sec)
    Loglikelihood: -20707210.979063
    Fst divergences between estimated populations: 
    	Pop0	Pop1	
    Pop0	
    Pop1	0.141	
    Pop2	0.120	0.128	
    Writing output files.


**Step 2**: Impute using ancestry estimates `P` and `Q`. Note `copyto!` and `convert` assumes `P` has dimension `K x S` and `Q` has dimension `K x N` where `K` is number of populations, `S` is number of SNPs, and `N` is number of individuals. So we need to transpose the output of `admixture`.


```julia
Pt = P |> transpose |> Matrix
Qt = Q |> transpose |> Matrix
convert(Matrix{Float64}, mouse, Pt, Qt)
```




    1940×10150 Matrix{Float64}:
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  …  2.0      2.0      2.0      2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     1.0      1.0      1.0      1.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  …  2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0     2.0      2.0      2.0      2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     1.0      1.0      1.0      1.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0  …  0.0      0.0      0.0      0.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     0.0      0.0      0.0      0.0
     ⋮                        ⋮         ⋱                             
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0  …  2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0     2.0      2.0      2.0      2.0
     2.0  2.0  2.0  2.0  2.0  2.0  2.0     2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  2.0  1.0  2.0  …  2.0      2.0      2.0      2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     1.0  1.0  2.0  1.0  1.0  1.0  1.0     2.0      2.0      2.0      2.0
     1.0  1.0  1.0  1.0  1.0  1.0  1.0     1.89365  1.89217  1.89208  1.89208
     0.0  0.0  0.0  0.0  2.0  0.0  2.0     2.0      2.0      2.0      2.0




```julia
# takes slightly longer because of calculation involving P and Q
M = similar(mouse, Float64)
@btime(copyto!($M, $mouse, $Pt, $Qt));
```

      138.702 ms (0 allocations: 0 bytes)


### Summaries

#### Counts

Counts of each the four possible values for each column are returned by `counts`.`


```julia
counts(mouse, dims=1)
```




    4×10150 Matrix{Int64}:
      358   359  252   358    33   359  …    56    56    56    56    56    56
        2     0    4     3     4     1      173   173   162   173   174   175
     1003  1004  888  1004   442  1004      242   242   242   242   242   242
      577   577  796   575  1461   576     1469  1469  1480  1469  1468  1467



Column 2 has no missing values (code `0x01`, the second row in the column-counts table).
In that SNP position for this sample, 359 indivduals are homozygous allele 1 (`G` according to the `.bim` file), 1004 are heterozygous, and 577 are homozygous allele 2 (`A`).

The counts by column and by row are cached in the `SnpArray` object. Accesses after the first are extremely fast.


```julia
@btime(counts($mouse, dims=1));
```

      4.698 ns (0 allocations: 0 bytes)


#### Minor allele frequencies

Minor allele frequencies (MAF) for each SNP.


```julia
maf(mouse)
```




    10150-element Vector{Float64}:
     0.4434984520123839
     0.4438144329896907
     0.359504132231405
     0.4439855446566856
     0.13119834710743805
     0.44404332129963897
     0.1412706611570248
     0.30299123259412064
     0.4445018069179143
     0.44424367578729995
     0.43427835051546393
     0.14075413223140498
     0.304639175257732
     ⋮
     0.0527624309392265
     0.052980132450331174
     0.08079096045197742
     0.08253250423968339
     0.08253250423968339
     0.10022650056625138
     0.10016977928692694
     0.10016977928692694
     0.09955005624296964
     0.10016977928692694
     0.10022650056625138
     0.10028328611898019



Minor allele (`false` means A1 is the minor allele; `true` means A2 is the minor allele) for each SNP.


```julia
minorallele(mouse)
```




    10150-element BitVector:
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     ⋮
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0



#### `mean` and `var`

The package provides methods for the generics `mean` and `var` from the `Statistics` package.


```julia
mean(mouse, dims=1)
```




    1×10150 Matrix{Float64}:
     1.113  1.11237  1.28099  1.11203  …  1.8009  1.79966  1.79955  1.79943




```julia
mean(mouse, dims=1, model=DOMINANT_MODEL)
```




    1×10150 Matrix{Float64}:
     0.815273  0.814948  0.869835  0.815178  …  0.968308  0.96829  0.968272




```julia
var(mouse, dims=1)
```




    1×10150 Matrix{Float64}:
     0.469929  0.470089  0.462605  0.469365  …  0.223714  0.223818  0.223923



These methods make use of the cached column or row counts and thus are very fast


```julia
@btime(mean($mouse, dims=1));
```

      15.272 μs (2 allocations: 79.39 KiB)


The column-wise or row-wise standard deviations are returned by `std`.


```julia
std(mouse, dims=2)
```




    1940×1 Matrix{Float64}:
     0.6504997290784408
     0.6379008244533891
     0.6558172726141286
     0.6532675479248437
     0.6744432174014563
     0.6519092298111158
     0.6779881845456428
     0.6955814098050999
     0.6437566832989493
     0.6505283141088536
     0.665444994623426
     0.659392039592328
     0.6641674726999468
     ⋮
     0.6599158250006595
     0.688387450736178
     0.6664063015924304
     0.6613451651895259
     0.6659810347614777
     0.6274577846909379
     0.6823658517777204
     0.6695299551061924
     0.710756592739754
     0.6387913736114869
     0.6736492722732016
     0.688855476425891



#### Missing rate

Proportion of missing genotypes


```julia
missingrate(mouse, 1)
```




    10150-element Vector{Float64}:
     0.0010309278350515464
     0.0
     0.002061855670103093
     0.0015463917525773195
     0.002061855670103093
     0.0005154639175257732
     0.002061855670103093
     0.0005154639175257732
     0.0015463917525773195
     0.0015463917525773195
     0.0
     0.002061855670103093
     0.0
     ⋮
     0.06701030927835051
     0.06597938144329897
     0.08762886597938144
     0.08814432989690722
     0.08814432989690722
     0.08969072164948454
     0.08917525773195877
     0.08917525773195877
     0.08350515463917525
     0.08917525773195877
     0.08969072164948454
     0.09020618556701031




```julia
missingrate(mouse, 2)
```




    1940-element Vector{Float64}:
     0.00019704433497536947
     0.0
     0.018423645320197045
     0.0007881773399014779
     0.0
     0.004236453201970443
     0.0051231527093596055
     0.00039408866995073894
     0.005517241379310344
     0.0016748768472906405
     0.0
     9.852216748768474e-5
     0.0004926108374384236
     ⋮
     0.000689655172413793
     0.004729064039408867
     0.0004926108374384236
     0.001083743842364532
     0.00019704433497536947
     0.0025615763546798028
     0.0038423645320197044
     0.001379310344827586
     0.0064039408866995075
     0.002857142857142857
     0.0011822660098522167
     0.00029556650246305416



#### Location of the missing values

The positions of the missing data are evaluated by


```julia
mp = missingpos(mouse)
```




    1940×10150 SparseArrays.SparseMatrixCSC{Bool, Int32} with 33922 stored entries:
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣯⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣽⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⠿⠿⠿⠿⠿⠿⠿⠺⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿⠿




```julia
@btime(missingpos($mouse));
```

      41.419 ms (19272 allocations: 1.80 MiB)


So, for example, the number of missing data values in each column can be evaluated as


```julia
sum(mp, dims=1)
```




    1×10150 Matrix{Int64}:
     2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175



although it is faster, but somewhat more obscure, to use


```julia
view(counts(mouse, dims=1), 2:2, :)
```




    1×10150 view(::Matrix{Int64}, 2:2, :) with eltype Int64:
     2  0  4  3  4  1  4  1  3  3  0  4  0  …  174  173  173  162  173  174  175



### Genetic relationship matrix (GRM)

#### Homogenous population

For homogenous population, `grm` function computes the empirical kinship matrix using either the classical genetic relationship matrix, `grm(A, model=:GRM)`, or the method of moment method, `grm(A, model=:MoM)`, or the robust method, `grm(A, model=:Robust)`. See the section _Kinship Comparison_ of the [manuscript](http://hua-zhou.github.io/media/pdf/Zhou19OpenMendel.pdf) for the formulae and references for these methods. 

Classical genetic relation matrix


```julia
# grm(mouse, method=:MoM)
# grm(mouse, method=:Robust)
g = grm(mouse, method=:GRM)
```




    1940×1940 Matrix{Float64}:
      0.478301    -0.0331304    0.0135612    …  -0.0347737   -0.0129443
     -0.0331304    0.422771    -0.0389227        0.0457987    0.00556832
      0.0135612   -0.0389227    0.509248        -0.0356689   -0.0608705
      0.0198205    0.00728645  -0.00935362      -0.0302404   -0.0102152
      0.056747    -0.0163418   -0.00495283      -0.0413347   -0.0415659
     -0.0165628   -0.0191127   -0.0112181    …   0.0177118   -0.0193087
      0.123771    -0.0404167    0.00442739       0.00880649  -0.0437565
     -0.0628362    0.172552    -0.0728312        0.0640027   -0.0281429
      0.0605018   -0.0260505    0.00398852      -0.00277754  -0.0607773
      0.108886    -0.0204594   -0.00767711      -0.0210501    0.00343526
     -0.0142307    0.00270989  -0.0235504    …  -0.0223563   -0.028408
     -0.0306022    0.197743    -0.00244269       0.0213998   -0.0478472
     -0.0131463   -0.0226707    0.0223522       -0.037288     0.0493662
      ⋮                                      ⋱               
      0.0176725   -0.0165609    0.0378308        0.0238751   -0.0420143
      0.0024949   -0.0411137    0.0154847       -0.0380656   -0.0650806
      0.0952286    0.00894298  -0.0163446    …  -0.0202633   -0.0219594
     -0.0309488   -0.0228342   -0.0478253       -0.014896     0.261623
     -0.004804    -0.0375168   -0.0211418       -0.0172572    0.0359166
      0.0076296    0.0481887   -0.0328968        0.0920425   -0.0292548
      0.070045    -0.0302138    0.000647283      0.00892069  -0.00632566
      0.0378132   -6.59565e-5   0.00888932   …   0.00230815  -0.0291622
     -0.00132837   0.00223654   0.0495928       -0.00936248   0.0299075
      0.0640864   -0.0241218    0.00602283       0.00403413   0.00689551
     -0.0347737    0.0457987   -0.0356689        0.509228    -0.035215
     -0.0129443    0.00556832  -0.0608705       -0.035215     0.552712




```julia
@btime(grm($mouse, method=:GRM));
```

      513.601 ms (15 allocations: 28.95 MiB)


Using Float32 (single precision) potentially saves memory usage and computation time.


```julia
grm(mouse, method=:GRM, t=Float32)
```




    1940×1940 Matrix{Float32}:
      0.478301    -0.0331304    0.0135612    …  -0.0347737   -0.0129443
     -0.0331304    0.422771    -0.0389227        0.0457987    0.00556833
      0.0135612   -0.0389227    0.509248        -0.0356689   -0.0608705
      0.0198205    0.00728645  -0.00935362      -0.0302404   -0.0102152
      0.056747    -0.0163418   -0.00495283      -0.0413347   -0.0415659
     -0.0165628   -0.0191127   -0.0112181    …   0.0177117   -0.0193087
      0.123771    -0.0404167    0.0044274        0.00880651  -0.0437565
     -0.0628362    0.172552    -0.0728312        0.0640027   -0.0281429
      0.0605018   -0.0260505    0.00398852      -0.00277754  -0.0607773
      0.108886    -0.0204594   -0.00767711      -0.0210501    0.00343525
     -0.0142307    0.0027099   -0.0235504    …  -0.0223563   -0.028408
     -0.0306022    0.197743    -0.00244268       0.0213998   -0.0478472
     -0.0131463   -0.0226707    0.0223522       -0.037288     0.0493662
      ⋮                                      ⋱               
      0.0176725   -0.016561     0.0378308        0.0238751   -0.0420143
      0.00249492  -0.0411137    0.0154847       -0.0380656   -0.0650806
      0.0952286    0.00894297  -0.0163446    …  -0.0202633   -0.0219594
     -0.0309488   -0.0228342   -0.0478253       -0.014896     0.261623
     -0.00480403  -0.0375167   -0.0211418       -0.0172572    0.0359166
      0.00762961   0.0481887   -0.0328968        0.0920425   -0.0292548
      0.070045    -0.0302138    0.000647267      0.0089207   -0.00632565
      0.0378132   -6.59479f-5   0.00888931   …   0.00230814  -0.0291622
     -0.00132835   0.00223653   0.0495928       -0.00936247   0.0299075
      0.0640864   -0.0241218    0.00602284       0.00403414   0.0068955
     -0.0347737    0.0457987   -0.0356689        0.509228    -0.035215
     -0.0129443    0.00556833  -0.0608705       -0.035215     0.552712




```julia
@btime(grm($mouse, method=:GRM, t=Float32));
```

      304.821 ms (16 allocations: 14.60 MiB)


By default, `grm` exlcudes SNPs with minor allele frequency below 0.01. This can be changed by the keyword argument `minmaf`.


```julia
# compute GRM excluding SNPs with MAF≤0.05 
grm(mouse, minmaf=0.05)
```




    1940×1940 Matrix{Float64}:
      0.478556    -0.0331783    0.013541     …  -0.0348225   -0.0129761
     -0.0331783    0.422993    -0.0389741        0.0457975    0.00554753
      0.013541    -0.0389741    0.50952         -0.0357183   -0.0609305
      0.0203209    0.00777944  -0.00887047      -0.0297696   -0.00972836
      0.0567523   -0.0163798   -0.00498406      -0.0413874   -0.0416146
     -0.0166009   -0.0191523   -0.0112531    …   0.0176939   -0.0193442
      0.123816    -0.0404689    0.00440171       0.0087834   -0.0438065
     -0.0629017    0.172626    -0.0729026        0.0640123   -0.0281836
      0.0605093   -0.0260942    0.00396257      -0.00280748  -0.0608373
      0.108922    -0.0204998   -0.00770996      -0.0210909    0.00341321
     -0.0142674    0.00268319  -0.0235927    …  -0.0223978   -0.0284489
     -0.0306486    0.197832    -0.00247243       0.0213842   -0.0478996
     -0.0131824   -0.0227124    0.0223371       -0.0373384    0.0493713
      ⋮                                      ⋱               
      0.0176546   -0.016599     0.0378249        0.0238609   -0.0420633
      0.00246808  -0.0411663    0.0154656       -0.0381165   -0.0651432
      0.0952566    0.00891997  -0.0163826    …  -0.0203036   -0.0219965
     -0.0309912   -0.0228718   -0.0478777       -0.0149289    0.261754
     -0.00483514  -0.0375673   -0.0211827       -0.0172957    0.0359138
      0.00770862   0.0482917   -0.0328417        0.0921714   -0.0292961
      0.0700582   -0.03026      0.000619365      0.00889767  -0.00635348
      0.0378313   -7.02155e-5   0.00889036   …   0.0023053   -0.0291795
     -0.00133338   0.00223364   0.0496179       -0.00937223   0.0299252
      0.0641201   -0.0241403    0.00602217       0.0040323    0.00689958
     -0.0348225    0.0457975   -0.0357183        0.509501    -0.0352599
     -0.0129761    0.00554753  -0.0609305       -0.0352599    0.553015



To specify specific SNPs for calculating empirical kinship, use the `cinds` keyword (default is `nothing`). When `cinds` is specified, `minmaf` is ignored.


```julia
# GRM using every other SNP
grm(mouse, cinds=1:2:size(mouse, 2))
```




    1940×1940 Matrix{Float64}:
      0.477       -0.0307774     0.0118026   …  -0.0320301    -0.0125113
     -0.0307774    0.425085     -0.0367459       0.0480442     0.00519065
      0.0118026   -0.0367459     0.505038       -0.0385129    -0.0631557
      0.0166017    0.00614789   -0.00919695     -0.0399744    -0.0104884
      0.05724     -0.0122148    -0.00543377     -0.0395663    -0.0372998
     -0.0193129   -0.0224378    -0.009277    …   0.0153785    -0.0220184
      0.12194     -0.0410682     0.00274307      0.00796748   -0.0441578
     -0.0624031    0.173985     -0.0724784       0.0663191    -0.0294243
      0.0627626   -0.0288615     0.00265615     -0.00449877   -0.0579702
      0.110878    -0.0232715    -0.00881604     -0.021272      0.00169016
     -0.00800735  -0.00149824   -0.019791    …  -0.024124     -0.0289397
     -0.0272944    0.19894      -0.00534771      0.0209384    -0.0511051
     -0.011388    -0.0281003     0.0273853      -0.0360047     0.0459359
      ⋮                                      ⋱                
      0.0169431   -0.0136989     0.0340794       0.0272811    -0.041189
      0.00201325  -0.0426611     0.0124353      -0.0387982    -0.0656181
      0.097587     0.0058123    -0.0160698   …  -0.021457     -0.023226
     -0.0342014   -0.0211246    -0.0490112      -0.0129575     0.256552
     -0.00324255  -0.0423482    -0.0192699      -0.0149015     0.0339388
      0.00575353   0.0464237    -0.0294694       0.0924759    -0.0275451
      0.0748725   -0.0258461    -0.00141068      0.0115232    -0.00486589
      0.0386555    0.000612169   0.00959997  …  -0.000357284  -0.0334687
     -0.00343056   0.0120673     0.0455375      -0.0103798     0.0336959
      0.0656909   -0.0193469     0.00600815      0.00188545    0.00726181
     -0.0320301    0.0480442    -0.0385129       0.513285     -0.0317963
     -0.0125113    0.00519065   -0.0631557      -0.0317963     0.54471



#### Inhomogenous/admixed populations

For inhomogenous/admixed population, we recommend first estimate the ancestry and pupulation allele frequencies using the ADMIXTURE software. See [ADMIXTURE.jl](https://github.com/OpenMendel/ADMIXTURE.jl) for usage. Then compute the kinship coefficients using the `P` (allele frequencies) and `Q` (ancestry fractions) matrix from the output of ADMIXTURE. This is essentially what the [REAP software](http://faculty.washington.edu/tathornt/software/REAP) does, except our implementation runs much faster than REAP (>50 fold speedup). 


```julia
# first read in the P and Q matrix output from ADMIXTURE and tranpose them
Pt = readdlm("mouse.3.P", ' ', Float64) |> transpose |> Matrix
Qt = readdlm("mouse.3.Q", ' ', Float64) |> transpose |> Matrix;
```


```julia
SnpArrays.grm_admixture(mouse, Pt, Qt)
```

    convert genotype: 0.26 seconds
    Φ = GG': 0.50 seconds
    convert G to {0,1} matrix: 0.02 seconds
    S = GG': 0.30 seconds





    1940×1940 Matrix{Float64}:
      0.459157     -0.0156933    -0.00323857  …  -0.0241032     0.0122942
     -0.0156933     0.38254      -0.00891405      0.00213638    0.0256826
     -0.00323857   -0.00891405    0.488813       -0.0171976    -0.0468329
      0.0030966     0.0202297    -0.0243853      -0.0176487     0.000561005
      0.0332113     0.00312454   -0.0272922      -0.032754     -0.0122725
     -0.0358108    -0.00853739   -0.0271515   …   0.0243009    -0.00138401
      0.121113     -0.0392918     0.00165846      0.010853     -0.0377001
     -0.0449199     0.10156      -0.0339266       0.0239008    -0.0152297
      0.0448067    -0.00136653   -0.0145566       0.00912369   -0.0448301
      0.0907814     0.0141711    -0.0255688      -0.00297642    0.00990265
      0.00337643    0.00227661   -0.00875794  …  -0.0258833    -0.0459813
     -0.0177263     0.145081      0.030058       -0.0133901    -0.0327255
      0.000957239  -0.00906892    0.029328       -0.0270566     0.0232713
      ⋮                                       ⋱                
     -0.00488414    0.000742334   0.0151452       0.0414038    -0.0193765
     -0.0206989    -0.0165344    -0.0109315      -0.0255689    -0.0508724
      0.0604748     0.0196283    -0.0426789   …  -0.0174937     0.027048
     -0.00510929   -0.00671179   -0.0330136      -0.00639737    0.21654
      0.0203133    -0.0184515    -0.0104571      -0.000598109  -0.0203991
      0.0143139     0.012726     -0.0201427       0.0734973    -0.0175892
      0.0431033    -0.0226698    -0.0277297       0.014573      0.0375869
      0.00723525    0.0106991    -0.0187392   …   0.00859215    0.00929037
     -0.0282601     0.0197691     0.0164641       0.0038119     0.0484681
      0.0376866    -0.00807382   -0.0201108       0.0180007     0.0385872
     -0.0241032     0.00213638   -0.0171976       0.501372     -0.0266089
      0.0122942     0.0256826    -0.0468329      -0.0266089     0.511513




```julia
# clean up
rm("mouse.3.P", force = true)
rm("mouse.3.Q", force = true)
```

### Filtering

Before GWAS, we often need to filter SNPs and/or samples according to genotyping success rates, minor allele frequencies, and Hardy-Weinberg Equilibrium test. This can be achieved by the `filter` function.

```@docs
SnpArrays.filter
```

By default, it outputs row and column index vectors such that sample-wise and SNP-wise genotyping success rate are at least 0.98 and minor allele frequencies are at least 0.01. User can opt to filter according to Hardy-Weinberg test by setting the minumum p-value `min_hwe_pval`.


```julia
rowmask, colmask =  SnpArrays.filter(mouse)
```




    (Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  1, 1, 1, 1, 1, 1, 1, 1, 1, 1], Bool[1, 1, 1, 1, 1, 1, 1, 1, 1, 1  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0])




```julia
count(rowmask), count(colmask)
```




    (1931, 10072)




```julia
@btime(SnpArrays.filter($mouse, min_success_rate_per_row=0.999, min_success_rate_per_col=0.999));
```

      135.726 ms (11459 allocations: 171.28 MiB)


One may use the `rowmask` and `colmask` to filter and save filtering result as Plink files.
```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask)
```

#### Filter Plink files

Filter a set of Plink files according to row indices and column indices. By result, filtered Plink files are saved as `srcname.filtered.bed`, `srcname.filtered.fam`, and `srcname.filtered.bim`, where `srcname` is the source Plink file name. You can also specify destimation file name using keyword `des`.


```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), 1:5, 1:5)
```




    5×5 SnpArray:
     0x02  0x02  0x02  0x02  0x03
     0x02  0x02  0x03  0x02  0x02
     0x03  0x03  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02
     0x03  0x03  0x03  0x03  0x03




```julia
# clean up
rm(SnpArrays.datadir("mouse.filtered.bed"), force=true)
rm(SnpArrays.datadir("mouse.filtered.fam"), force=true)
rm(SnpArrays.datadir("mouse.filtered.bim"), force=true)
```

Filter a set of Plink files according to logical vectors.


```julia
SnpArrays.filter(SnpArrays.datadir("mouse"), rowmask, colmask)
```




    1931×10072 SnpArray:
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x02  0x03  0x02  0x02  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x03  0x00  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x00  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x03  0x02  0x02  0x02  0x03
     0x02  0x02  0x02  0x02  0x03  0x02  …  0x03  0x00  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x00  0x03  0x00  0x00
     0x02  0x02  0x03  0x02  0x02  0x02     0x00  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x00  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x03  0x02  0x02  0x02  0x03
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x00  0x03  0x00  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x00  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x00  0x00  0x03  0x03  0x03
        ⋮                             ⋮  ⋱                             ⋮  
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x00  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x03  0x00  0x00  0x00
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x00  0x03  0x00  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x00  0x03
     0x03  0x03  0x03  0x03  0x03  0x03  …  0x02  0x03  0x03  0x00  0x03  0x03
     0x02  0x02  0x02  0x02  0x03  0x02     0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x00  0x00  0x03  0x00  0x03
     0x02  0x02  0x03  0x02  0x02  0x02     0x03  0x03  0x00  0x03  0x03  0x03
     0x02  0x02  0x02  0x02  0x02  0x02     0x03  0x00  0x00  0x03  0x03  0x03
     0x00  0x00  0x00  0x00  0x03  0x00  …  0x03  0x03  0x00  0x03  0x00  0x03




```julia
readdir(glob"mouse.filtered.*", datapath)
```




    3-element Vector{String}:
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.filtered.bed"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.filtered.bim"
     "/home/xyz/.julia/dev/SnpArrays/data/mouse.filtered.fam"




```julia
# clean up
rm(SnpArrays.datadir("mouse.filtered.bed"), force=true)
rm(SnpArrays.datadir("mouse.filtered.fam"), force=true)
rm(SnpArrays.datadir("mouse.filtered.bim"), force=true)
```

### Concatenating `SnpArray`s

Concatenation of `SnpArray`s is implemented in `hcat`, `vcat`, and `hvcat` functions. By default, the resulting `.bed` file is saved as a file beginning with `tmp_` in the working directory. You can specify destination using keyword `des`. 

For concatenation, `SnpArray` arguments do not deal with `.fam` or `.bim` files at all. You can use `SnpData` as the arguments to create those files (see below).


```julia
s = SnpArrays.filter(SnpArrays.datadir("mouse"), 1:2, 1:3)
s
```




    2×3 SnpArray:
     0x02  0x02  0x02
     0x02  0x02  0x03




```julia
all(s .== [[0x02 0x02 0x02];
[0x02 0x02 0x03]])
```




    true



Standard concatenation works just like any other arrays. However, a temporary file is created as a side effect.


```julia
[s s s]
```




    2×9 SnpArray:
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03




```julia
[s; s; s]
```




    6×3 SnpArray:
     0x02  0x02  0x02
     0x02  0x02  0x03
     0x02  0x02  0x02
     0x02  0x02  0x03
     0x02  0x02  0x02
     0x02  0x02  0x03




```julia
[s s s; s s s]
```




    4×9 SnpArray:
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03




```julia
readdir(glob"tmp_*", ".")
```




    7-element Vector{String}:
     "./tmp_hcat_arr_1.bed"
     "./tmp_hvcat_arr_1.bed"
     "./tmp_vcat_arr_1.bed"
     "./tmp_vcat_arr_2.bed"
     "./tmp_vcat_arr_3.bed"
     "./tmp_vcat_arr_4.bed"
     "./tmp_vcat_arr_5.bed"



In order to set the destination `.bed` file, you can add the keyword argument `des`.


```julia
hcat(s, s, s; des=SnpArrays.datadir("mouse.test.hcat"))
```




    2×9 SnpArray:
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03




```julia
vcat(s, s, s; des=SnpArrays.datadir("mouse.test.vcat"))
```




    6×3 SnpArray:
     0x02  0x02  0x02
     0x02  0x02  0x03
     0x02  0x02  0x02
     0x02  0x02  0x03
     0x02  0x02  0x02
     0x02  0x02  0x03




```julia
hvcat((3, 3), s, s, s, s, s, s; des=SnpArrays.datadir("mouse.test.hvcat"))
```




    4×9 SnpArray:
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03
     0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02  0x02
     0x02  0x02  0x03  0x02  0x02  0x03  0x02  0x02  0x03




```julia
# clean up
rm(SnpArrays.datadir("mouse.filtered.bed"), force=true)
rm(SnpArrays.datadir("mouse.filtered.fam"), force=true)
rm(SnpArrays.datadir("mouse.filtered.bim"), force=true)
tmplist = readdir(glob"tmp_*.bed", ".")
for f in tmplist
    rm(f, force=true)
end
rm(SnpArrays.datadir("mouse.test.hcat.bed"), force=true)
rm(SnpArrays.datadir("mouse.test.vcat.bed"), force=true)
rm(SnpArrays.datadir("mouse.test.hvcat.bed"), force=true)
```

## Linear Algebra

In some applications we want to perform linear algebra using SnpArray directly without expanding it to numeric matrix. This is achieved in three different `struct`s:

1. Direct operations on a plink-formatted `SnpArray`: `SnpLinAlg`
2. Operations on transformed `BitMatrix`es: `SnpBitMatrix`
3. Direct operations on a plink-formatted data on an Nvidia GPU: `CuSnpArray`.

`SnpLinAlg` and `SnpBitMatrix` use Chris Elrod's [LoopVectorization.jl](https://github.com/chriselrod/LoopVectorization.jl) internally. It is much faster on machines with AVX support. `CuSnpArray` uses [CUDA.jl](https://juliagpu.gitlab.io/CUDA.jl/) internally.

!!! warning "deprecated SnpBitMatrix"
    `SnpBitMatrix` is now deprecated in favor of `SnpLinAlg`. 
    `SnpBitMatrix` will be removed on next minor release.
    
The implementation assumes that the matrix corresponding to SnpArray is the matrix of the A2 allele counts. `SnpLinAlg` and `CuSnpArray` impute any missing genotype with its column mean by default. They can also configured to impute missing genotypes with zero. `SnpBitMatrix` can only impute missing values with zero. 

### Constructor

First let's load a data set without missing genotypes.


```julia
const EUR = SnpArray(SnpArrays.datadir("EUR_subset.bed"))
```




    379×54051 SnpArray:
     0x03  0x03  0x03  0x02  0x02  0x03  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x02  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x03  0x03  0x03     0x02  0x02  0x02  0x03  0x03  0x02
     0x03  0x03  0x03  0x00  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x03  0x03     0x02  0x02  0x02  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x03  0x03  0x03  0x02
     0x02  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x02  0x02  0x03  0x03
     0x02  0x03  0x03  0x03  0x02  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x03  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03  …  0x03  0x03  0x02  0x02  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x02
     0x03  0x02  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
        ⋮                             ⋮  ⋱     ⋮                             ⋮
     0x03  0x03  0x03  0x00  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x02  0x02  0x02  0x03  0x02  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x03  0x03  …  0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x00  0x00  0x03     0x02  0x02  0x02  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x03  0x03  0x03  …  0x03  0x03  0x02  0x02  0x03  0x03
     0x03  0x03  0x03  0x00  0x03  0x03     0x03  0x03  0x03  0x03  0x03  0x03
     0x02  0x03  0x03  0x02  0x00  0x02     0x03  0x03  0x03  0x03  0x03  0x03
     0x03  0x03  0x03  0x02  0x02  0x03     0x03  0x03  0x03  0x03  0x03  0x03



To instantiate a SnpLinAlg based on SnpArray,


```julia
const EURsla = SnpLinAlg{Float64}(EUR, model=ADDITIVE_MODEL, center=true, scale=true)
const EURsla_ = SnpLinAlg{Float64}(EUR, model=ADDITIVE_MODEL, center=true, scale=true, impute=false)
```




    379×54051 SnpLinAlg{Float64}:
      0.46516   0.163517  0.306468  -0.0298581  …   0.342518   0.163517   0.23281
      0.46516  -6.0338    0.306468  -0.0298581      0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468   1.38467        0.342518   0.163517  -4.17894
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468   1.38467    …   0.342518   0.163517  -4.17894
     -1.91722   0.163517  0.306468  -0.0298581     -2.7483     0.163517   0.23281
     -1.91722   0.163517  0.306468   1.38467        0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581  …  -2.7483     0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518   0.163517  -4.17894
      0.46516  -6.0338    0.306468  -0.0298581      0.342518   0.163517   0.23281
      ⋮                                         ⋱                         ⋮
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518  -6.0338     0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468  -0.0298581  …   0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468   1.38467        0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468   1.38467    …  -2.7483     0.163517   0.23281
      0.46516   0.163517  0.306468  -1.44439        0.342518   0.163517   0.23281
     -1.91722   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281
      0.46516   0.163517  0.306468  -0.0298581      0.342518   0.163517   0.23281



The constructor shares the same keyword arguments as the `convert` or `copyto!` functions. The type parameter, `Float64` in this example, indicates the SnpLinAlg acts like a Float64 matrix.
SnpLinAlg directly uses the SnpArray for computation and does not expand into full numeric array. 


```julia
Base.summarysize(EUR), Base.summarysize(EURsla)
```




    (6876757, 8609701)



### `mul!`

SnpLinAlg act similar to a regular matrix and responds to `size`, `eltype`, SnpLinAlg-vector multiplication, and SnpLinAlg-matrix multiplications. Other linear algebra operations (e.g. `qr()`) should work on a SnpLinAlg, but will be much slower. 


```julia
@show size(EURsla)
@show eltype(EURsla)
@show typeof(EURsla) <: AbstractMatrix;
```

    size(EURsla) = (379, 54051)
    eltype(EURsla) = Float64
    typeof(EURsla) <: AbstractMatrix = true


Matrix-vector and matrix-matrix multiplications with SnpLinAlg are mathematically equivalent to the corresponding Float matrix contained from `convert` or `copyto!` a SnpArray.


```julia
using LinearAlgebra
v1 = randn(size(EUR, 1))
v2 = randn(size(EUR, 2))
A = convert(Matrix{Float64}, EUR, model=ADDITIVE_MODEL, center=true, scale=true);
```


```julia
norm(EURsla * v2 - A * v2)
```




    3.4786140310420274e-11




```julia
norm(EURsla' * v1 - A' * v1)
```




    5.471334812385415e-12



### Linear Algebra Performance

See Linear Algebra Benchmarks on the left for performance comparison among BLAS, SnpLinAlg, and CuSnpArray (for GPU). In general,
+ SnpLinAlg-vector multiplications are at least 2x faster than the corresponding Matrix{Float64}-vector multiplication using BLAS
+ CuSnpArray-vector multiplications on the GPU is 50x faster than BLAS, and
+ SnpLinAlg-matrix multiplication is competitive with BLAS if the right hand matrix is "tall and thin".

Note that SnpLinAlg does not allocate additional memory, and can impute missing values with column means. 

### `copyto!`, `convert`, and `subarrays`

`copyto!` and `convert` are also supported on `SnpLinAlg`s, but without the `impute`, `scale`, `center` keyword arguments. The destination array will be scaled/centered if the `SnpLinAlg` was scaled/centered. 


```julia
# convert work on SnpLinAlg (and subarrays of it)
Atrue = convert(Matrix{Float64}, EUR, center=true, scale=true, impute=true)
A = convert(Matrix{Float64}, EURsla)
all(Atrue .≈ A)
```




    true




```julia
# copyto on a subarray
v = zeros(size(EUR, 1), 10)
copyto!(v, @view(EURsla[:, 1:2:20]))
```




    379×10 Matrix{Float64}:
      0.46516  0.306468  -0.539104  -0.370294  …   0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
     -1.91722  0.306468   0.97438   -0.370294  …   0.238721   0.551318   0.301719
     -1.91722  0.306468  -0.539104  -1.83218      -4.06963    0.551318   0.301719
     -1.91722  0.306468  -0.539104  -0.370294      0.238721  -1.53818    0.301719
      0.46516  0.306468  -0.539104  -0.370294     -4.06963    0.551318   0.301719
     -1.91722  0.306468   0.97438   -0.370294      0.238721   0.551318  -3.16348
      0.46516  0.306468   0.97438    1.09159   …   0.238721  -1.53818    0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
      0.46516  0.306468  -0.539104   1.09159       0.238721   0.551318   0.301719
      ⋮                                        ⋱                        
      0.46516  0.306468  -0.539104  -0.370294      0.238721   0.551318   0.301719
      0.46516  0.306468  -0.539104  -0.370294      0.238721   0.551318   0.301719
      0.46516  0.306468  -0.539104  -0.370294     -4.06963    0.551318   0.301719
     -1.91722  0.306468   0.97438   -0.370294  …   0.238721   0.551318   0.301719
      0.46516  0.306468  -2.05259   -0.370294      0.238721   0.551318   0.301719
     -1.91722  0.306468   0.97438   -0.370294      0.238721  -1.53818    0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
     -1.91722  0.306468   0.97438   -0.370294  …   0.238721   0.551318   0.301719
      0.46516  0.306468   0.97438    1.09159       0.238721   0.551318   0.301719
     -1.91722  0.306468  -2.05259   -1.83218       0.238721   0.551318   0.301719
      0.46516  0.306468  -0.539104  -0.370294      0.238721   0.551318   0.301719




```julia
all(v .≈ Atrue[:, 1:2:20])
```




    true



### GPU support: CuSnpArray (optional)

On machines with Nvidia GPU, matrix-vector multiplications can be performed on it via CuSnpArray. The input vectors should be CuVectors. 


```julia
using CUDA, Adapt
out1 = randn(size(EUR, 1))
out2 = randn(size(EUR, 2))
v1 = randn(size(EUR, 1))
v2 = randn(size(EUR, 2))
v1_d = adapt(CuVector{Float64}, v1) # sends data to GPU
v2_d = adapt(CuVector{Float64}, v2)
out1_d = adapt(CuVector{Float64}, out1)
out2_d = adapt(CuVector{Float64}, out2)

const EURcu = CuSnpArray{Float64}(EUR; model=ADDITIVE_MODEL, center=true, scale=true);
```

    ┌ Warning: The NVIDIA driver on this system only supports up to CUDA 10.2.0.
    │ For performance reasons, it is recommended to upgrade to a driver that supports CUDA 11.2 or higher.
    └ @ CUDA /home/xyz/.julia/packages/CUDA/CtvPY/src/initialization.jl:42



```julia
@btime mul!($out1_d, $EURcu, $v2_d);
```


```julia
@btime mul!($out2_d, transpose($EURcu), $v1_d);
```

The operations are parallelized along the output dimension, hence the GPU was not fully utilized in the first case. With 100-time larger data, 30 to 50-fold speedup were observed for both cases with Nvidia Titan V. See linear algebra page for more information.

Let's check correctness of the result.

## SnpData

We can create a `SnpData`, which has a `SnpArray` with information on SNP and subject appended.

### Constructor


```julia
EUR_data = SnpData(SnpArrays.datadir("EUR_subset"))
```

### Filter

We can filter SnpData by functions `f_person` and `f_snp`. `f_person` applies to the field `person_info` and selects persons (rows) for which `f_person` is `true`.`f_snp` applies to the field `snp_info` and selects snps (columns) for which `f_snp` is `true`. The first argument can be either a `SnpData` or an `AbstractString`.


```julia
SnpArrays.filter(EUR_data; des="tmp.filter.chr.17", f_snp = x -> x[:chromosome]=="17")
```


```julia
SnpArrays.filter(SnpArrays.datadir("EUR_subset"); des="tmp.filter.chr.17", f_snp = x -> x[:chromosome]=="17")
```


```julia
SnpArrays.filter(EUR_data; des="tmp.filter.sex.male", f_person = x -> x[:sex] == "1")
```

Both `f_person` and `f_snp` can be used at the same time.


```julia
SnpArrays.filter(EUR_data; des="tmp.filter.chr.17.sex.male", f_person = x -> x[:sex] == "1", f_snp = x -> x[:chromosome] == "17")
```

### Split

We can split `SnpData` by SNP's choromosomes or each person's sex or phenotype using `split_plink`. Again, the first argument can be an `SnpData` or an `AbstractString`.


```julia
splitted = SnpArrays.split_plink(SnpArrays.datadir("EUR_subset"), :chromosome; prefix="tmp.split.chr.")
```

Let's take a SnpArray for chromosome 17.


```julia
piece = splitted["17"]
```


```julia
@assert all(piece.snp_info[!, :chromosome].== "17")
```


```julia
splitted_sex = SnpArrays.split_plink(EUR_data, :sex; prefix="tmp.split.sex.")
```

### Concatenation

`hcat`, `vcat`, and `hvcat` are also implemented for `SnpData`. All of `.bed`, `.bim`, `.fam` files are created. Simple concatenation expression can be used (with the side effect of creation of temporary plink files). One may also set the desitination using the keyword argument `des`. 


```julia
[piece piece]
```


```julia
[piece; piece]
```


```julia
[piece piece; piece piece]
```


```julia
hcat(piece, piece; des="tmp.hcat")
```


```julia
vcat(piece, piece; des="tmp.vcat")
```


```julia
hvcat((2,2), piece, piece, piece, piece; des="tmp.hvcat")
```

### Merge

We can merge the splitted dictionary back into one SnpData using `merge_plink`.


```julia
merged = SnpArrays.merge_plink("tmp.merged", splitted) # write_plink is included here
```

You can also merge the plink formatted files based on their common prefix.


```julia
merged_from_splitted_files = merge_plink("tmp.split.chr"; des = "tmp.merged.2")
```

### Reorder

Order of subjects can be changed using the function `reorder!`.


```julia
const mouse_prefix = SnpArrays.datadir("mouse")
run(`cp $(mouse_prefix * ".bed") mouse_reorder.bed`)
run(`cp $(mouse_prefix * ".bim") mouse_reorder.bim`)
run(`cp $(mouse_prefix * ".fam") mouse_reorder.fam`)
```


```julia
mouse_data = SnpData(mouse_prefix)
mouse_toreorder = SnpData("mouse_reorder", "r+")
m, n = size(mouse_toreorder.snparray)
```

For example, the below randomly permutes subjects.


```julia
using Random
ind = randperm(m)
SnpArrays.reorder!(mouse_toreorder, ind)
```


```julia
mouse_toreorder
```

This functionality mainly targets Cox regression, where sorting subjects in decreasing order of (censored) survival time results in more efficient implementation.

## VCF to PLINK

SnpArrays.jl includes a function to transform a (gzipped) VCF file to PLINK-formatted files. This function drops multi-allelic variants and variants with missing identifier.


```julia
# Download an example VCF file
isfile("test.08Jun17.d8b.vcf.gz") || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz", 
    joinpath(pwd(), "test.08Jun17.d8b.vcf.gz"));
```


```julia
vcf2plink("test.08Jun17.d8b.vcf.gz", "test.08Jun17.d8b")
```


```julia
# clean up
for ft in ["bim", "fam", "bed"]
    rm("tmp.filter.chr.17." * ft, force=true)
    rm("tmp.filter.sex.male." * ft, force=true)
    rm("tmp.filter.chr.17.sex.male." * ft, force=true)
    for k in keys(splitted)
        rm("tmp.split.chr.$(k)." * ft, force=true)
    end
    for k in keys(splitted_sex)
        rm("tmp.split.sex.$(k)." * ft, force=true)
    end
    rm("tmp.merged." * ft, force=true)
    rm("tmp.merged.2." * ft, force=true)
    
    rm("tmp.hcat." * ft, force=true)
    rm("tmp.vcat." * ft, force=true)
    rm("tmp.hvcat." * ft, force=true)

    tmplist = glob("tmp_*" * ft)
    for f in tmplist
        rm(f, force=true)
    end
end
tmplist = readdir(glob"tmp_*.bed", ".")
for f in tmplist
    rm(f, force=true)
end
rm("mouse_reorder.bim", force=true)
rm("mouse_reorder.bed", force=true)
rm("mouse_reorder.fam", force=true)
rm("mouse_reorder.reordered.fam", force=true)
rm("test.08Jun17.d8b.vcf.gz", force=true)
rm("test.08Jun17.d8b.bed", force=true)
rm("test.08Jun17.d8b.bim", force=true)
rm("test.08Jun17.d8b.fam", force=true)
```
